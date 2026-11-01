import numpy as np


class Cell(np.ndarray):
    """
    Array containing any kind of object (inherits from ndarray).
    The data type is defined as `numpy.object`.
    This is similar to cells in Matlab.
    """

    def __new__(subtype, shape, buffer=None, offset=0,
                strides=None, order=None):
        dtype = object
        obj = super(Cell, subtype).__new__(subtype, shape, dtype,
                                           buffer, offset, strides,
                                           order)
        return obj


def isCell(u):
    """
    Test if object `u` is an instance of Cell class.

    Parameters
    ----------
    u : any object
            
    Returns
    -------
    bool
        True if 'u' is an instance of Cell class, False otherwise.
    """
    return isinstance(u, Cell)
