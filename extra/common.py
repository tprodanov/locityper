import sys
import gzip
import numpy as np

def create_power_mean(power):
    if power == 1:
        return np.mean
    elif power == np.inf:
        return np.max
    elif power == -np.inf:
        return np.min

    import scipy
    return functools.partial(scipy.stats.pmean, p=power)


def open_stream(filename, mode='r'):
    assert mode == 'r' or mode == 'w'
    if filename is None or filename == '-':
        return sys.stdin if mode == 'r' else sys.stdout
    elif filename.endswith('.gz'):
        return gzip.open(filename, mode + 't')
    else:
        return open(filename, mode)
