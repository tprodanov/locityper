import sys
import gzip
import lzma
import os
import builtins
import csv
import subprocess
try:
    from compression.zstd import open as zstd_open
except ModuleNotFoundError:
    def zstd_open(*args, **kwargs):
        error('zstd compression not supported (required Python >= 3.14)')


class BrotliReader:
    def __init__(self, filename, mode='rt'):
        assert 'w' not in mode
        self.proc = subprocess.Popen(
            ['brotli', '-Kdc', filename],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text='b' not in mode,
        )

    def close(self):
        proc = self.proc
        self.proc = None
        proc.stdout.close()
        stderr = proc.stderr.read()
        proc.stderr.close()
        ret = proc.wait()
        if ret:
            raise subprocess.CalledProcessError(ret, proc.args, stderr=stderr)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()

    def __getattr__(self, name):
        return getattr(self.proc.stdout, name)


class BrotliWriter:
    def __init__(self, filename, mode='wt'):
        assert 'r' not in mode
        self.proc = subprocess.Popen(
            ['brotli', '-5', '-f', '-w', '24', '-o', filename],
            stdin=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text='b' not in mode,
        )

    def close(self):
        proc = self.proc
        self.proc = None
        proc.stdin.close()
        stderr = proc.stderr.read()
        proc.stderr.close()
        ret = proc.wait()
        if ret:
            raise subprocess.CalledProcessError(ret, proc.args, stderr=stderr)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()

    def __getattr__(self, name):
        return getattr(self.proc.stdin, name)


def open(filename, mode='rt'):
    if filename is None or filename == '-':
        return sys.stdin if 'r' in mode else sys.stdout
    elif filename.endswith('.gz'):
        return gzip.open(filename, mode)
    elif filename.endswith('.xz'):
        return lzma.open(filename, mode)
    elif filename.endswith('.zst'):
        level = 7 if 'w' in mode else None
        return zstd_open(filename, mode, level=level)
    elif filename.endswith('.br'):
        return BrotliReader(filename, mode) if 'r' in mode else BrotliWriter(filename, mode)
    else:
        return builtins.open(filename, mode)


def mkdir(path):
    try:
        os.mkdir(path)
    except FileExistsError:
        pass


def read_csv(f):
    """
    Wrapper over csv.DictReader, that ignores first lines with comments.
    """
    fields = None
    for line in f:
        if not line.startswith('#'):
            fields = line.strip().split('\t')
            break
    assert fields is not None
    return csv.DictReader(f, fieldnames=fields, delimiter='\t')


def error(msg):
    sys.stderr.write(f'ERROR: {msg}\n')
    exit(1)


class Sink:
    def write(self, text):
        pass
