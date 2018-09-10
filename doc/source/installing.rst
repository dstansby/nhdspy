Installing
==========

To install nhdspy, download the code from https://github.com/dstansby/nhdspy
Then change into the source directory and run::

  pip install nhdspy

To run nhdspy you will need the following requirements:

  - gfortran
  - hdf5
  - libz
  - szip

Running nhdspy will automatically try to compile the bundled FORTRAN code.
The following environment variables must point to their respective libaries:

  - HDF5
  - LIBZ
  - LIBSZ
