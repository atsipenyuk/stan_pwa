import numpy as np

def array_to_string(a):
    """
    Convert the numpy array to string.

    The returned string may be later executed as python code to define 
    the numpy array.
    """
    nr = a.shape[0]
    nc = a.shape[1]

    res = "np.asarray(["
    for i in range(nr):

        res += "["
        for k in range(nc):
            res += str(a[i,k]) + ","

        # Delete the last comma and close the row
        res = res[:-1] + "],"

    # Delete the last comma and close the array
    res = res[:-1] + "])"
    return res



def vector_to_string(a):
    """Convert numpy vector to string."""

    nr = a.shape[0]

    res = "np.asarray(["
    for i in range(nr):
        res += str(a[i]) + ","

    # Delete the last comma and close the array
    res = res[:-1] + "])"
    return res
