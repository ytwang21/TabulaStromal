"""Atomic artifacts, content signatures, validated checkpoints and process locks."""
from pathlib import Path
from contextlib import contextmanager
import hashlib, json, os, tempfile, time
import joblib
import numpy as np
import pandas as pd
_HASHES = {}


def json_default(x):
    if isinstance(x, Path): return str(x)
    if isinstance(x, np.ndarray): return x.tolist()
    if isinstance(x, np.generic): return x.item()
    if isinstance(x, set): return sorted(x)
    raise TypeError(type(x).__name__)


def digest(x):
    return hashlib.sha256(json.dumps(x, sort_keys=True, default=json_default).encode()).hexdigest()


def file_info(path):
    p = Path(path)
    return dict(path=str(p.resolve()), size=p.stat().st_size, mtime_ns=p.stat().st_mtime_ns)


def file_hash(path):
    p = Path(path)
    st = p.stat()
    # Atomic replacement can preserve size and timestamp ticks on a coarse
    # filesystem. A new inode must never reuse the prior file's digest.
    key = (str(p.resolve()), st.st_dev, st.st_ino, st.st_size, st.st_mtime_ns, st.st_ctime_ns)
    if key in _HASHES: return _HASHES[key]
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1024 * 1024), b''): h.update(b)
    result = h.hexdigest()
    _HASHES[key] = result
    return result


@contextmanager
def atomic_path(path):
    p = Path(path); p.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix='.' + p.name + '.', suffix='.tmp', dir=p.parent)
    os.close(fd)
    try:
        yield Path(tmp)
        with open(tmp, 'rb') as f: os.fsync(f.fileno())
        os.replace(tmp, p)
    finally:
        if os.path.exists(tmp): os.unlink(tmp)


def write_json(path, obj):
    with atomic_path(path) as tmp:
        tmp.write_text(json.dumps(obj, indent=2, default=json_default, allow_nan=True) + '\n')


def write_text(path, text):
    with atomic_path(path) as tmp: tmp.write_text(text)


def write_csv(path, df):
    with atomic_path(path) as tmp: df.to_csv(tmp, index=False)


def dump(path, obj):
    with atomic_path(path) as tmp: joblib.dump(obj, tmp, compress=3)


class Cache:
    def __init__(self, signature, resume=True):
        self.signature, self.resume = signature, resume
        self.created = set()

    def valid(self, marker):
        p = Path(marker)
        if (not self.resume and str(p.resolve()) not in self.created) or not p.exists(): return False
        try:
            d = json.loads(p.read_text())
            if d.get('signature') != self.signature or d.get('status') != 'complete': return False
            for f in d.get('outputs', []):
                child = Path(f['path'])
                if not child.is_file() or child.stat().st_size != f['size'] or file_hash(child) != f['sha256']: return False
                if child.name.endswith('complete.json') and not self.valid(child): return False
            return True
        except (OSError, ValueError, KeyError):
            return False

    def mark(self, marker, outputs, **extra):
        write_json(marker, dict(signature=self.signature, status='complete', completed_at=time.time(),
                               outputs=[dict(path=str(Path(p).resolve()), size=Path(p).stat().st_size,
                                             sha256=file_hash(p)) for p in outputs], **extra))
        self.created.add(str(Path(marker).resolve()))


@contextmanager
def process_lock(path):
    """fcntl releases the lock even after SIGKILL; lock files need not be deleted."""
    import fcntl
    p = Path(path); p.parent.mkdir(parents=True, exist_ok=True)
    with p.open('a+') as f:
        try: fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError: raise RuntimeError(f'Another pipeline owns {p}')
        f.seek(0); f.truncate(); f.write(str(os.getpid())); f.flush()
        yield
