#This file is part of NHDS.
#
#    NHDS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    NHDS is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with NHDS.  If not, see <http://www.gnu.org/licenses/>.

SRCS := $(shell ls src/*.f90)
OBJS := $(SRCS:src/%.f90=obj/%.o)

HDF5 := /usr/local/Cellar/hdf5/1.10.2_1
LIBZ := /usr/local/Cellar/zlib/1.2.11
LIBSZ := /usr/local/Cellar/szip/2.1.1_1

INCLUDEHDF := -I$(HDF5)/include
LIBSHDF := $(HDF5)/lib/libhdf5_fortran.a $(HDF5)/lib/libhdf5.a
LIBSLIBZ := $(LIBZ)/lib/libz.a
LIBSLIBSZ := $(LIBSZ)/lib/libsz.a

COMPILER := gfortran -fbounds-check -O2 -ffixed-line-length-none -ffree-line-length-none -Wunused -lm $(INCLUDEHDF)


bin/NHDS: $(OBJS)
	$(COMPILER) $(OBJS) -o $@ $(LIB64) $(LIBSHDF) $(LIBSLIBZ) $(LIBSLIBSZ)


obj/%.o:src/%.f90
	$(COMPILER)  -c  $< -o $@


clean:
	rm -f $(OBJS) bin/NHDS
