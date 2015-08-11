###############################################################################
#                                                                             #
# Convert python dict to STAN *.data.R functions.                             #
#                                                                             #
# The functions below make a STAN file from python dict.                      #
# Copied from STAN Github project with some minor simplifications.            #
#                                                                             #
# SOURCE: https://github.com/stan-dev/pystan/blob/develop/pystan/misc.py.     #
#                                                                             #
###############################################################################

import numpy as np
import os
import ROOT

def _dict_to_rdump(data):
    parts = []
    for name, value in data.items():
        #if isinstance(value, (np.number, np.ndarray, int, bool, float)) \
        #    and not isinstance(value, string_types):
        value = np.asarray(value)
        #else:
        #    raise ValueError("Variable {} is not a number and cannot be dumped.".format(name))

        #if value.dtype == np.bool:
        #    value = value.astype(int)

        if value.ndim == 0:
            s = '{0} <- {1}\n'.format(name, str(value))
        elif value.ndim == 1:
            if len(value) == 1:
                s = '{0} <-\nc({1})\n'.format(name, ', '.join(str(v) for v in value))#                s = '{0} <- {1}\n'.format(name, str(value.pop()))
            else:
                s = '{0} <-\nc({1})\n'.format(name, ', '.join(str(v) for v in value))
        elif value.ndim > 1:
            tmpl = '{0} <-\nstructure(c({1}), .Dim = c({2}))\n'
            s = tmpl.format(name,
                            ', '.join(str(v) for v in value.flatten(order='F')),
                            ', '.join(str(v) for v in value.shape))
        parts.append(s)
    return ''.join(parts)


def stan_rdump(data, filename):
    """
    Dump a dictionary with model data into a file using the R dump format that
    Stan supports.

    Parameters
    ----------
    data : dict
    filename : str
    """
    #for name in data:
    #    if not is_legal_stan_vname(name):
    #        raise ValueError("Variable name {} is not allowed in Stan".format(name))
    with open(filename, 'w') as f:
        f.write(_dict_to_rdump(data))

# END OF GITHUB STAN CODE #####################################################
