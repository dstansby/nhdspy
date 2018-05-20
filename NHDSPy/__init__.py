import subprocess
import tempfile
import numpy as np
from . import helpers


class Species:
    """
    A single bi-Maxwellian species.

    Parameters
    ----------
    q : float
        Particle charge, in units of proton charge. e.g. for a proton
        enter ``1``, for an electron enter ``-1``.
    m : float
        Particle mass, in units of proton mass.
    n : float
        Number density, relative to the number density used to define the
        Alfvén speed.
    v_d : float
        Drift speed
    t_ani : float
        Temperature anisotropy (perpendicular temperature divided by
        parallel temperature)
    beta_par : float
        Parallel beta.
    """
    def __init__(self, q, m, n, v_d, t_ani, beta_par):
        assert m > 0, 'Mass must be greater than zero'
        self.q = q
        self.m = m
        self.n = n
        self.v_d = v_d
        self.t_ani = t_ani
        self.beta_par = beta_par


class InputParams:
    """
    Input parameters for calculating a dispersion relation.

    Parameters
    ----------
    species : list of :class:`Species`
        List of the particle species for which to calculate the dispersion.
    propagation_angle : float
        Propagation angle to calculate dispersion relation for. Defined
        with respect to the magnetic field.
    va : float
        Ratio of the Alfvén speed to the speed of light.

    Other parameters
    ----------------
    numiter : int, optional
        Maximum number of iterations in the Newton method
    det_D_threshold : float, optional
        Threshold for the determinant of the dispersion tensor.
        If det D <= det_D_threshold, the Newton iteration will be stopped.
    n_bessel : int, optional
        Maximum of sum in Bessel functions (both regular and modified).
        Can be very low (e.g., 3) for quasi-parallel propagation.
    bessel_zero : float, optional
        If I_n is less than this value, higher n are neglected.
    vxsteps, vysteps, vzsteps : int, optional
        Steps in the x,y,z directions.
    """
    def __init__(self, species, propagation_angle, va,
                 numiter=1000, det_D_threshold=1e-16,
                 n_bessel=1000, bessel_zero=1e-50,
                 vxsteps=100, vysteps=100, vzsteps=100
                 ):
        self.species = species
        self.propagation_angle = propagation_angle
        self.va = va
        self.numiter = numiter
        self.det_D_threshold = det_D_threshold
        self.n_bessel = n_bessel
        self.bessel_zero = bessel_zero
        self.vxsteps = vxsteps
        self.vysteps = vysteps
        self.vzsteps = vzsteps

    @property
    def nspecies(self):
        '''
        Number of species.
        '''
        return len(self.species)


class Result:
    def __init__(self, input, run_output):
        self.input = input
        self.output = self._process_output(run_output)

    def _process_output(self, run_output):
        split = run_output.split('\n')[:-1]
        split = [[float(a) for a in s.split()] for s in split]
        return np.array(split)

    @property
    def kz(self):
        return self.output[:, 0]

    @property
    def omega_real(self):
        '''
        Real part of frequency normalised to the proton gyro-frequency.
        '''
        return self.output[:, 1]

    @property
    def omega_imag(self):
        '''
        Imaginary part of frequency normalised to the proton gyro-frequency.
        '''
        return self.output[:, 2]


def format_input_file(input):
    '''
    Function to create input file.

    Parameters
    ----------
    input : InputParams
    '''
    input_file_sting = r'''
    subroutine set_parameters()
    use input_params
    implicit none

    ! Number of species
    ! (up to 10 species are possible)
    numspec={numspec}

    ! Maximum number of iterations in the Newton method
    numiter={numiter}

    ! Threshold for the determinant of the dispersion tensor:
    ! If det D <= det_D_threshold, the Newton iteration will be stopped
    det_D_threshold={detD}

    ! Maximum of sum in Bessel functions (both regular and modified)
    ! can be very low (e.g., 3) for quasi-parallel propagation
    nmax={nmax}

    ! If I_n is less than this value, higher n are neglected:
    Bessel_zero={bessel_zero}

    ! Temperature anisotropy (Tperp/Tparallel)
    alpha=(/ 1.d0,1.d0,1.d0,1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0 /)

    ! Parallel beta of the species
    beta=(/ 1.d0,1.d0,1.d0,1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0 /)

    ! Charge of the species in units of the first ion charge
    charge=(/ 1.d0,-1.d0,2.d0,-1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0 /)

    ! Mass of the species in units of ion mass
    mass=(/ 1.d0,1.d0/1836.d0,4.d0,1.d0/1836.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0 /)

    ! Density of the species in units of proton density
    density=(/ 1.d0,1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0 /)

    ! Drift speed of the species in units of proton Alfven speed
    vdrift=(/ 0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0 /)

    ! Angle of propagation (in degrees)
    theta={propagation_angle}

    ! Alfven speed divided by speed of light
    vAc={va}

    ! Amplitude mode: If 1, then ampl = delta B_x / B0. If 2, then ampl = delta B_y / B0.
    ! If 3, then ampl = delta B_z / B0
    ampl_mode=3

    ! Amplitude for the calculation of polarization properties:
    ampl=1.d0

    ! Print warnings to std. output and stop program in case of problems
    output_warning=.FALSE.

    !! The following parameters are used to determine delta f in write_delta_f:
    ! Maximum of sum in Bessel function:
    ! can be very low (e.g., 3) for quasi-parallel propagation
    mmax=1000

    ! If I_n is less than this value, higher m are neglected:
    Bessel_zero_deltaf=1.d-50

    ! Steps in vx, vy, and vz (standard: 100):
    vxsteps={vxsteps}
    vysteps={vysteps}
    vzsteps={vzsteps}

    ! Range in vpar and vperp:
    vxrange=(/ -1.d0,1.d0 /)
    vyrange=(/ -1.d0,1.d0 /)
    vzrange=(/ -1.d0,1.d0 /)

    ! Number of time steps for one full period (standard: 25):
    timesteps=40

    ! If periods is .TRUE., then num_periods is number of periods (2 PI / omega_r).
    ! If periods is .FALSE., then num_periods is number of gyro-periods (2 PI / Omega_p).
    periods=.TRUE.

    ! Number of periods (standard: 8):
    num_periods=1

    ! Include wave damping in delta f calculation or not (exp(gamma*t)-term):
    damping=.FALSE.

    ! If const_r is .TRUE., then delta f is evaluated at point r = 0 -> cos(-omega * t)
    ! If const_r is .FALSE., then the wave period is a function of space at -> cos(k * r)
    !		with k * r from 0 ... 2 pi at a fixed time. We "fly" along the k-direction over one wave train.
    ! 		If set on .FALSE., periods does not make a difference
    const_r=.TRUE.

    end subroutine
    '''.format(numspec=input.nspecies,
               numiter=input.numiter,
               detD=helpers._float_to_fortran_str(input.det_D_threshold),
               nmax=input.n_bessel,
               bessel_zero=helpers._float_to_fortran_str(input.bessel_zero),
               propagation_angle=helpers._float_to_fortran_str(input.propagation_angle),
               va=helpers._float_to_fortran_str(input.va),
               vxsteps=input.vxsteps,
               vysteps=input.vysteps,
               vzsteps=input.vzsteps,
               )
    return input_file_sting


def run(input):
    '''
    Run the dispersion solver for a given input.

    Parameters
    ----------
    input : InputParams
    '''
    NHDS_binary = 'NHDS/bin/NHDS'
    out = subprocess.check_output([NHDS_binary], universal_newlines=True)
    out = Result(input, out)
    print(out.output.shape)
    return out
