import sys
import gzip
import os
import builtins


def open(filename, mode='r'):
    assert mode == 'r' or mode == 'w'
    if filename is None or filename == '-':
        return sys.stdin if mode == 'r' else sys.stdout
    elif filename.endswith('.gz'):
        return gzip.open(filename, mode + 't')
    else:
        return builtins.open(filename, mode)


def mkdir(path):
    try:
        os.mkdir(path)
    except FileExistsError:
        pass
