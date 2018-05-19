import numpy as np


def _float_to_fortran_str(fl):
    if fl == 0:
        exponent = 0
        factor = 0.0
    else:
        exponent = int(np.floor(np.log10(fl)))
        factor = fl / 10**exponent
    return '{}d{}'.format(factor, exponent)
