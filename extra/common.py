import sys
import gzip
import lzma
import os
import builtins
import csv
from contextlib import contextmanager
import subprocess


@contextmanager
def read_brotli(filename, mode='rt'):
    proc = subprocess.Popen(
        ['brotli', '-Kdc', filename],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text='b' not in mode,
    )
    try:
        yield proc.stdout
    finally:
        proc.stdout.close()
        stderr = proc.stderr.read()
        proc.stderr.close()
        ret = proc.wait()
        if ret:
            raise subprocess.CalledProcessError(ret, proc.args, stderr=stderr)


@contextmanager
def write_brotli(filename, mode='wt'):
    proc = subprocess.Popen(
        ['brotli', '-5', '-f', '-w', '24', '-o', filename],
        stdin=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text='b' not in mode,
    )
    try:
        yield proc.stdin
    finally:
        proc.stdin.close()
        stderr = proc.stderr.read()
        proc.stderr.close()
        ret = proc.wait()
        if ret:
            raise subprocess.CalledProcessError(ret, proc.args, stderr=stderr)


def open(filename, mode='rt'):
    if filename is None or filename == '-':
        return sys.stdin if 'r' in mode else sys.stdout
    elif filename.endswith('.gz'):
        return gzip.open(filename, mode)
    elif filename.endswith('.xz'):
        return lzma.open(filename, mode)
    elif filename.endswith('.br'):
        return read_brotli(filename, mode) if 'r' in mode else write_brotli(filename, mode)
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
