import subprocess
import tempfile
import numpy as np


class Species:
    """

    """
    def __init__(self, q, m, n, v_d, t_ani, beta_par):
        pass


class Input:
    """
    Input parameters for calculating a dispersion relation.

    Parameters
    ----------
    numiter : int, optional
        Maximum number of iterations in the Newton method
    det_D_threshold : float, optional
        Threshold for the determinant of the dispersion tensor.
        If det D <= det_D_threshold, the Newton iteration will be stopped.
    bessel_zero : float, optional
        If I_n is less than this value, higher n are neglected.
    """
    def __init__(self, numiter=1000, det_D_threshold=1e-16, Bessel_zero=1e-50):
        self.numiter = numiter
        self.det_D_threshold = det_D_threshold


class Result:
    def __init__(self, run_output):
        self.output = self._process_output(run_output)

    def _process_output(self, run_output):
        split = run_output.split('\n')[:-1]
        split = [[float(a) for a in s.split()] for s in split]
        return np.array(split)


def run():
    NHDS_binary = 'NHDS/bin/NHDS'
    out = subprocess.check_output([NHDS_binary], universal_newlines=True)
    out = Result(out)
    print(out.output.shape)
