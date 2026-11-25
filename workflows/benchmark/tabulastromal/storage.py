"""Atomic artifacts, content signatures, validated checkpoints and process locks."""
from pathlib import Path
from contextlib import contextmanager
import hashlib, json, os, tempfile, time
import joblib
import numpy as np
import pandas as pd
_HASHES = {}
_WORKER_LOG = None


def worker_logging(root,phase):
    import logging,sys
    global _WORKER_LOG
    key=(os.getpid(),str(root),phase)
    if _WORKER_LOG==key:return
    _WORKER_LOG=key
    logging.basicConfig(level=logging.INFO,force=True,
        format='%(asctime)s | %(process)d | %(levelname)s | %(name)s | %(message)s',
        handlers=[logging.FileHandler(Path(root)/'logs'/f'{phase}_{os.getpid()}.log'),logging.StreamHandler(sys.stdout)])


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
    return dict(path=str(p.resolve()), size=p.stat().st_size, mtime_ns=p.stat().st_mtime_ns,ctime_ns=p.stat().st_ctime_ns)


def annotation_identity(path, cfg):
    """Hash actual annotation/category codes and gene IDs without loading expression.

    Stat identity also covers matrix replacement. Annotation hashes explicitly
    invalidate old pseudobulk/model checkpoints after Celltype relabeling.
    """
    import h5py
    h=hashlib.sha256()
    with h5py.File(path,'r') as f:
        fields=set(v for v in cfg['fields'].values() if isinstance(v,str))
        fields.update(cfg['fields'].get('age_fallbacks',[]));fields.update(['_index','tissue','sex'])
        nodes=[f'obs/{x}' for x in sorted(fields) if x in f['obs']]
        nodes += [x for x in ['var/_index','raw/var/_index'] if x in f]
        def add(node):
            h.update(node.name.encode())
            if isinstance(node,h5py.Group):
                for key in sorted(node):add(node[key])
            else:
                for start in range(0,node.shape[0],65536):
                    arr=node[start:start+65536]
                    if arr.dtype.kind in 'OSU':h.update(json.dumps([str(v) for v in arr]).encode())
                    else:h.update(arr.tobytes())
        for node in nodes:add(f[node])
    return h.hexdigest()


def file_hash(path):
    p = Path(path)
    st = p.stat()
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
        self.started = time.time()

    def valid(self, marker):
        p = Path(marker)
        if not p.exists(): return False
        try:
            d = json.loads(p.read_text())
            # A worker's successful marker belongs to this fresh run even when
            # its process-local `created` set cannot be returned to the driver.
            if not self.resume and str(p.resolve()) not in self.created and d.get('completed_at',0)<self.started:return False
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
def directory_lock(path):
    """Lock the existing preprocessing directory: no extra lock files."""
    import fcntl
    p=Path(path);p.mkdir(parents=True,exist_ok=True);fd=os.open(p,os.O_RDONLY)
    try:
        fcntl.flock(fd,fcntl.LOCK_EX)
        yield
    finally:
        fcntl.flock(fd,fcntl.LOCK_UN);os.close(fd)


@contextmanager
def process_lock(path,blocking=False):
    """fcntl releases the lock even after SIGKILL; lock files need not be deleted."""
    import fcntl
    p = Path(path); p.parent.mkdir(parents=True, exist_ok=True)
    with p.open('a+') as f:
        try: fcntl.flock(f, fcntl.LOCK_EX | (0 if blocking else fcntl.LOCK_NB))
        except BlockingIOError: raise RuntimeError(f'Another pipeline owns {p}')
        f.seek(0); f.truncate(); f.write(str(os.getpid())); f.flush()
        yield
