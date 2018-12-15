# This file is part of NHDS
# Copyright (C) 2018 Daniel Verscharen (d.verscharen@ucl.ac.uk)
#All rights reserved.
#
#Redistribution and use in source and binary forms, with or without
#modification, are permitted provided that the following conditions are met:
#
#1. Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#2. Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
#ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
#(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
#LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
#ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#The views and conclusions contained in the software and documentation are those
#of the authors and should not be interpreted as representing official policies,
#either expressed or implied, of the NHDS project.

SRCS := $(shell ls src/*.f90)
OBJS := $(SRCS:src/%.f90=obj/%.o)

# INCLUDEHDF := -I$(HDF5)/include
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
