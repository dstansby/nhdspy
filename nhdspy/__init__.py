import os
import pathlib as path
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
        Particle mass, in units of proton mass. e.g. for a proton
        enter ``1``.
    n : float
        Number density, relative to the number density used to define the
        Alfvén speed.
    v_d : float
        Drift speed as a fraction of the Alfvén speed.
    t_ani : float
        Temperature anisotropy (perpendicular temperature divided by
        parallel temperature).
    beta_par : float
        Parallel beta (thermal pressure divided by magnetic pressure).
    """
    def __init__(self, q, m, n, v_d, t_ani, beta_par):
        assert m > 0, 'Mass must be greater than zero'
        self.q = q
        self.m = m
        self.n = n
        self.v_d = v_d
        self.t_ani = t_ani
        self.beta_par = beta_par

    def __repr__(self):
        return ('Species object with ' +
                'q={self.q}, m={self.m}, n={self.n}, '.format(self=self) +
                'v_d={self.v_d}, t_ani={self.t_ani}, beta_par={self.beta_par}'.format(self=self)
                )


class InputParams:
    """
    Input parameters for calculating a dispersion relation.

    Parameters
    ----------
    species : list of :class:`Species`
        List of the particle species for which to calculate the dispersion.
        **Note** that the order of the species matters, as some of the
        parameters below are defined relative to the first species in the
        *species* list.
    omega_guess : complex
        Initial guess for the frequency, normalised to the gyro frequency of
        the first species.
    propagation_angle : float
        Propagation angle to calculate dispersion relation for. Defined
        with respect to the magnetic field.
    va : float
        Ratio of the Alfvén speed to the speed of light. The number density
        in the Alfvén speed is taken from the first species.
    kzmin : float
        Start of kz range.
    kzmax : float
        End of kz range.

    Other parameters
    ----------------
    kzsteps : int, optional
        Number of steps in kz to calculate dispersion relation at. Points are
        linearly spaced between *kzmin* and *kzmax*.
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
    def __init__(self, species, omega_guess, propagation_angle, va,
                 kzmin, kzmax, kzsteps=200,
                 numiter=1000, det_D_threshold=1e-16,
                 n_bessel=1000, bessel_zero=1e-50,
                 vxsteps=100, vysteps=100, vzsteps=100
                 ):
        self.species = species
        self.omega_guess = omega_guess
        self.propagation_angle = propagation_angle
        self.va = va
        self.kzmin = kzmin
        self.kzmax = kzmax
        self.kzsteps = kzsteps
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

    @property
    def total_charge(self):
        """
        Total charge density. Result is normalised to the first species' charge
        and number density.
        """
        q = 0
        for species in self.species:
            q += species.q * species.n
        return q

    @property
    def total_current(self):
        """
        Total current density.Result is normalised to the first species'
        charge, number density, and the Alfvén speed.
        """
        j = 0
        for species in self.species:
            j += species.q * species.n * species.v_d
        return j


class Result:
    '''
    Result of running the dispersion solver.
    '''
    def __init__(self, input, run_output):
        self.input = input
        self.output = run_output

    @property
    def kz(self):
        '''
        Wave vector normalised to the proton intertial length (:math:`kd_{p}`).
        '''
        return self.output[:, 0]

    @property
    def omega_real(self):
        r'''
        Real part of frequency normalised to the proton gyro-frequency
        (:math:`\omega / \Omega_{p}`).
        '''
        return self.output[:, 2]

    @property
    def omega_imag(self):
        r'''
        Imaginary part of frequency normalised to the proton gyro-frequency
        (:math:`\gamma / \Omega_{p}`).
        '''
        return self.output[:, 3]

    @property
    def EyEx(self):
        '''
        The electric field component ratio :math:`E_{y} / E_{x}`.
        Returned as complex numbers.
        '''
        return self.output[:, 4] + self.output[:, 5] * 1j

    @property
    def EzEx(self):
        '''
        The electric field component ratio :math:`E_{z} / E_{x}`.
        Returned as complex numbers.
        '''
        return self.output[:, 6] + self.output[:, 7] * 1j


def format_input_file(input):
    '''
    Function to create input file.

    Parameters
    ----------
    input : InputParams
    '''
    def create_species_list(attr):
        '''
        Function to populate species properties.
        '''
        n = 0
        ret = ''
        for species in input.species:
            x = getattr(species, attr)
            ret += '{},'.format(helpers._float_to_fortran_str(x))
            n += 1
        for i in range(n, 10):
            ret += '0.d0,'
        ret = ret[:-1]
        return ret

    input_file_sting = r'''
&parameters

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

! Initial guess for the frequency (complex number):
initial_guess=({omega_guess_real},{omega_guess_imag})

! Range of values to scan over in kz:
krange={kzmin},{kzmax}

! Number of steps to scan over kzrange:
ksteps={kzsteps}

! Temperature anisotropy (Tperp/Tparallel)
alpha={alpha}

! Parallel beta of the species
beta={beta}

! Charge of the species in units of the first ion charge
charge={charge}

! Mass of the species in units of ion mass
mass={mass}

! Density of the species in units of proton density
density={density}

! Drift speed of the species in units of proton Alfven speed
vdrift={vdrift}

! Range of angle of propagation (in degrees)
theta_range={propagation_angle}, {propagation_angle}

! Steps in the theta direction
theta_steps=1

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
vxrange=-1.d0,1.d0
vyrange=-1.d0,1.d0
vzrange=-1.d0,1.d0

! Number of time steps for one full period (standard: 25):
timesteps=40

! If periods is .TRUE., then num_periods is number of periods (2 PI / omega_r).
! If periods is .FALSE., then num_periods is number of gyro-periods (2 PI / Omega_p).
periods=T

! Number of periods (standard: 8):
num_periods=1

! Include wave damping in delta f calculation or not (exp(gamma*t)-term):
damping=F

! If const_r is .TRUE., then delta f is evaluated at point r = 0 -> cos(-omega * t)
! If const_r is .FALSE., then the wave period is a function of space at -> cos(k * r)
!		with k * r from 0 ... 2 pi at a fixed time. We "fly" along the k-direction over one wave train.
! 		If set on .FALSE., periods does not make a difference
const_r=T

! Write moments to a file or not:
output_mom=F

! Write E & B vectors to a file or not:
output_EB=F

! Override the k,theta values specified above with the ones read from a file:
kth_file=F

! Filename to read k, theta values from:
! kth_filename='kth_kawtc'
/
    '''.format(numspec=input.nspecies,
               numiter=input.numiter,
               detD=helpers._float_to_fortran_str(input.det_D_threshold),
               nmax=input.n_bessel,
               bessel_zero=helpers._float_to_fortran_str(input.bessel_zero),
               omega_guess_real=helpers._float_to_fortran_str(input.omega_guess.real),
               omega_guess_imag=helpers._float_to_fortran_str(input.omega_guess.imag),
               propagation_angle=helpers._float_to_fortran_str(input.propagation_angle),
               va=helpers._float_to_fortran_str(input.va),
               vxsteps=input.vxsteps,
               vysteps=input.vysteps,
               vzsteps=input.vzsteps,
               alpha=create_species_list('t_ani'),
               beta=create_species_list('beta_par'),
               charge=create_species_list('q'),
               mass=create_species_list('m'),
               density=create_species_list('n'),
               vdrift=create_species_list('v_d'),
               kzmin=helpers._float_to_fortran_str(input.kzmin),
               kzmax=helpers._float_to_fortran_str(input.kzmax),
               kzsteps=input.kzsteps
               )
    return input_file_sting


def run(input):
    '''
    Run the dispersion solver for a given input.

    Parameters
    ----------
    input : InputParams
    '''
    abspath = path.Path(os.path.dirname(__file__))
    nhds_folder = abspath / 'NHDS'
    (nhds_folder / 'obj').mkdir(exist_ok=True)
    (nhds_folder / 'bin').mkdir(exist_ok=True)
    parameter_file = nhds_folder / 'parameters.in'
    # Create input file
    with open(parameter_file, mode='w') as f:
        f.write(format_input_file(input))
    # Check if NHDS is compiled
    binary = nhds_folder / 'NHDS'
    if not binary.exists():
        helpers.compile_nhds(nhds_folder)
    subprocess.check_output(
        ['./NHDS', 'parameters.in'], universal_newlines=True, cwd=nhds_folder)

    output_file = nhds_folder / 'output_parameters.in.dat'
    out = np.loadtxt(output_file, skiprows=1)
    return Result(input, out)
