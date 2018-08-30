Installing
==========

To install NHDSPy, download the code from https://github.com/dstansby/NHDSPy
Then change into the source directory and run::

  pip install NHDSPy

To run NHDSPy you will need the following requirements:

  - gfortran
  - hdf5
  - libz
  - szip

Running NHDSPy will automatically try to compile the bundled FORTRAN code.
The following environment variables must point to their respective libaries:

  - HDF5
  - LIBZ
  - LIBSZ
