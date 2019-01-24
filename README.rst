nhdspy
======

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.2302057.svg
   :target: https://doi.org/10.5281/zenodo.2302057

.. image:: https://readthedocs.org/projects/nhdspy/badge/?version=latest
   :target: https://nhdspy.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

.. image:: https://travis-ci.org/dstansby/nhdspy.svg?branch=master
   :target: https://travis-ci.org/dstansby/nhdspy

.. image:: https://codecov.io/gh/dstansby/nhdspy/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/dstansby/nhdspy


nhdspy is a wrapper for the NHDS dispersion relation solver. For the original NHDS code see https://github.com/danielver02/NHDS.

If you use this software to produce data for publication, please cite the NHDS paper: http://iopscience.iop.org/article/10.3847/2515-5172/aabfe3

Installing
----------

To install nhdspy, run

.. code-block:: bash

  pip install nhdspy


Example
-------

.. code-block:: python

 import nhdspy
 import matplotlib.pyplot as plt

 electrons = nhdspy.Species(-1, 1 / 1836, 1, 0, 1, 1)
 protons = nhdspy.Species(1, 1, 1, 0, 1, 1)
 propagation_angle = 0.001

 input = nhdspy.InputParams([protons, electrons], propagation_angle, 1e-4)
 print(nhdspy.format_input_file(input))

 output = nhdspy.run(input)

 fig, ax = plt.subplots()
 ax.plot(output.kz, output.omega_real)
 ax.plot(output.kz, output.omega_imag)
 plt.show()
