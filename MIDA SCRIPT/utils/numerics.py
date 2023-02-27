"""
Numerical utilities.

`binnings` solves the item binning problem. `cartesian` is a fast cartesian
product.

Author: Casey W. Stark <caseywstark@gmail.com>
Affiliation: UC Berkeley
Homepage: http://caseywstark.com
License:
  Copyright (C) 2011, 2012 Casey W. Stark. All Rights Reserved.

  This file is part of `MIDA`.

"""

import copy
import itertools

import numpy as np

# from http://stackoverflow.com/questions/6750298/efficient-item-binning-algorithm-itertools-numpy
def binnings(items, bins, cache={}):
    """
    Generate an array of all item binning possibilities.

    Parameters
    ----------
    items : int
        Number of items to put in bins.
    bins : int
        Number of bins.
    cache : dictionary
        The cache for recursion.

    Returns
    -------
    result : ndarray
        2-D array of shape (C(bins + items - 1, items), bins). I think...

    """
    # catch possible ends of recursion
    if items == 0:
        return np.zeros((1, bins), dtype=np.int32)
    if bins == 0:
        return np.empty((0, 0), dtype=np.int32)

    # setup key and check cache for this one
    args = (items, bins)
    if args in cache:
        return cache[args]

    # black magic goes here
    a = binnings(items - 1, bins, cache)

    # add 1 to 0-th column
    a1 = a + (np.arange(bins) == 0)
    b = binnings(items, bins - 1, cache)

    # insert zeros 0-th column
    b1 = np.hstack((np.zeros((b.shape[0], 1), dtype=np.int32), b))
    result = np.vstack((a1, b1))

    # don't forget to cache this one
    cache[args] = result

    return result

def binnings_iterator(items, bins):
    """
    Returns an iterator of all item binning possibilities. Note that the rows
    are just tuples, not ndarrays.

    """
    return itertools.ifilter(lambda combo: sum(combo) == items,
                             itertools.product(xrange(items + 1), repeat=bins))

# from http://stackoverflow.com/questions/1208118/using-numpy-to-build-an-array-of-all-combinations-of-two-arrays
def cartesian(arrays, out=None):
    """
    Generate a cartesian product of input arrays.

    Parameters
    ----------
    arrays : list of array-like
        1-D arrays to form the cartesian product of.
    out : numpy.ndarray
        Array to place the cartesian product in.

    Returns
    -------
    out : ndarray
        2-D array of shape (M, len(arrays)) containing cartesian products
        formed of input arrays.

    Examples
    --------
    >>> cartesian(([1, 2, 3], [4, 5], [6, 7]))
    array([[1, 4, 6],
           [1, 4, 7],
           [1, 5, 6],
           [1, 5, 7],
           [2, 4, 6],
           [2, 4, 7],
           [2, 5, 6],
           [2, 5, 7],
           [3, 4, 6],
           [3, 4, 7],
           [3, 5, 6],
           [3, 5, 7]])

    Note that the function fails if the arrays are too large.

    """
    arrays = [np.asarray(x) for x in arrays]
    dtype = arrays[0].dtype

    n = np.prod([x.size for x in arrays])
    if out is None:
        try:
            out = np.zeros([n, len(arrays)], dtype=dtype)
        except ValueError:
            # this happens when the array is too big.
            # fall back on iterator method.
            print "HIT big array."

            # @todo: try to iterate it out and make an array
            return
            #limited_prod = itertools.ifilter(lambda c: (c * np.arange(c.shape[1])).sum() <= 50,
            #                                 itertools.product(arrays))
            #return np.array(list(limited_prod), dtype=np.int)

    m = n / arrays[0].size
    out[:,0] = np.repeat(arrays[0], m)
    if arrays[1:]:
        cartesian(arrays[1:], out=out[0:m,1:])
        for j in xrange(1, arrays[0].size):
            out[j*m:(j+1)*m,1:] = out[0:m,1:]
    return out
