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


def temporary_filename(filename):
    return os.path.join(os.path.dirname(filename), 'tmp-{}'.format(os.path.basename(filename)))
