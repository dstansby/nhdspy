import os
import subprocess

import numpy as np


def _float_to_fortran_str(fl):
    if fl == 0:
        exponent = 0
        factor = 0.0
    else:
        exponent = int(np.floor(np.log10(np.abs(fl))))
        factor = fl / 10**exponent
    return '{}d{}'.format(factor, exponent)


def compile_nhds(nhds_folder):
    subprocess.run(['make', 'clean'], cwd=nhds_folder)
    result = subprocess.run(['make'], cwd=nhds_folder)
    if result.returncode != 0:
        print('')
        raise RuntimeError(
            'Compiling NHDS failed, see above for error messages')
    os.rename((nhds_folder / 'bin' / 'NHDS'), nhds_folder / 'NHDS')
