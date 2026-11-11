"""Read plotting metadata from the supplied protocol-4 AnnData pickles.

The original expression matrices are not required for these figures. Large
NumPy buffers are therefore skipped by seeking over their on-disk payloads.
The resulting AnnData objects are deliberately metadata-only and must never
be used for expression analysis or saved as complete AnnData files.
"""
from dataclasses import dataclass
from pathlib import Path
import io
import pickle
import struct
import time

import numpy as np


@dataclass
class DeferredBytes:
    offset: int
    length: int
    inline: bytes | None = None


class MetadataUnpickler(pickle._Unpickler):
    """Plot-only loader; ordinary metadata and UMAP arrays remain unchanged."""
    dispatch = pickle._Unpickler.dispatch.copy()

    def __init__(self, stream, threshold=16*1024**2):
        super().__init__(stream)
        self.stream = stream
        self.threshold = threshold
        self.deferred_arrays = {}
        self.skipped_bytes = 0
        self.start_time = time.monotonic()
        self.last_message = self.start_time

    def buffer(self, length):
        if length < self.threshold:
            self.append(self.read(length))
            return
        frame = self._unframer.current_frame
        if frame is not None:
            remaining = len(frame.getbuffer()) - frame.tell()
            if remaining == 0:
                self._unframer.current_frame = None
            elif remaining >= length:
                self.append(DeferredBytes(-1, length, frame.read(length)))
                return
            else:
                raise pickle.UnpicklingError('Large buffer crosses a pickle frame')
        offset = self.stream.tell()
        self.stream.seek(length, io.SEEK_CUR)
        self.append(DeferredBytes(offset, length))
        self.skipped_bytes += length
        now = time.monotonic()
        if now-self.last_message >= 15:
            print(f'Pickle metadata: reached {self.stream.tell()/1024**3:.1f} GiB; '
                  f'skipped {self.skipped_bytes/1024**3:.1f} GiB expression/unused buffers', flush=True)
            self.last_message = now

    def load_binbytes(self):
        self.buffer(struct.unpack('<I', self.read(4))[0])

    def load_binbytes8(self):
        self.buffer(struct.unpack('<Q', self.read(8))[0])

    def payload(self, buffer):
        if buffer.inline is not None:
            return buffer.inline
        position = self.stream.tell()
        self.stream.seek(buffer.offset)
        raw = self.stream.read(buffer.length)
        self.stream.seek(position)
        if len(raw) != buffer.length:
            raise pickle.UnpicklingError('Truncated deferred array')
        return raw

    def load_build(self):
        state, instance = self.stack[-1], self.stack[-2]
        if isinstance(instance, np.ndarray) and isinstance(state, tuple) and len(state)==5 and isinstance(state[-1], DeferredBytes):
            self.stack.pop()
            version, shape, dtype, fortran, buffer = state
            # Preserve UMAP (n x 2), and possible large pandas metadata blocks
            # (few fields x n observations). Expression arrays stay deferred.
            if len(shape)==2 and (shape[-1]==2 or shape[0]<=128):
                instance.__setstate__((version, shape, dtype, fortran, self.payload(buffer)))
            else:
                self.deferred_arrays[id(instance)] = (instance, shape, dtype, fortran, buffer)
            return
        super().load_build()

    dispatch[pickle.BINBYTES[0]] = load_binbytes
    dispatch[pickle.BINBYTES8[0]] = load_binbytes8
    dispatch[pickle.BUILD[0]] = load_build

    def array(self, value):
        """Materialize only a selected plotting array, never an expression X."""
        spec = self.deferred_arrays.get(id(value))
        if spec is None:
            return np.asarray(value)
        _, shape, dtype, fortran, buffer = spec
        raw = self.payload(buffer)
        return np.frombuffer(raw, dtype=dtype).reshape(shape, order='F' if fortran else 'C').copy()


def normalize_records(obj):
    """Confirm actual pickle structure here when using another input source."""
    import pandas as pd
    if isinstance(obj, pd.DataFrame):
        return obj.to_dict('records')
    if isinstance(obj, (list, tuple)):
        return list(obj)
    if isinstance(obj, dict):
        if 'adata' in obj or 'recall_df' in obj:
            return [obj]
        result=[]
        for key, value in obj.items():
            if not isinstance(value, dict):
                raise TypeError('Dataset dictionary values must be record dictionaries')
            result.append(dict(value, _dataset_key=key))
        return result
    raise TypeError(f'Unsupported top-level pickle type: {type(obj)}')


def metadata_records(path, full=False):
    """Yield raw records and a loader while the source file remains open."""
    from contextlib import contextmanager

    @contextmanager
    def opened():
        with Path(path).open('rb') as stream:
            if full:
                obj = pickle.load(stream)
                yield normalize_records(obj), None
            else:
                loader=MetadataUnpickler(stream)
                obj=loader.load()
                print(f'Read metadata from {Path(path).name}: '
                      f'{loader.skipped_bytes/1024**3:.2f} GiB large buffers skipped', flush=True)
                yield normalize_records(obj), loader
    return opened()
