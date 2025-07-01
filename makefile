# SPDX-License-Identifier: MIT
# This file is licensed under the MIT License.
# See https://opensource.org/licenses/MIT

FC   = gfortran
# FOPT = -ffree-line-length-none -fopenmp -fdefault-real-8 -z execstack -J build/ -O0 -g3 -ggdb -fcheck=all -fbacktrace -finit-integer=-9999 -finit-real=-inf
FOPT = -ffree-line-length-none -fopenmp -fdefault-real-8 -z execstack -J build/ -O2

CC = gcc
COPT = -z execstack -O2
CLIB = -lgfortran -lm -lquadmath -lgomp

all: build/example build/distribution_table.o build/generalized_SG.a build/C_interface.a build/C_example

build/util.o: util.f90
	${FC} ${FOPT} -c util.f90 -o build/util.o

build/quad.o: quad.f90
	${FC} ${FOPT} -c quad.f90 -o build/quad.o

build/distribution_table.o: distribution_table.f90 build/quad.o build/util.o
	${FC} ${FOPT} -c distribution_table.f90 -o build/distribution_table.o

build/generalized_SG.o: generalized_SG.f90 build/distribution_table.o build/quad.o build/util.o
	${FC} ${FOPT} -c generalized_SG.f90 -o build/generalized_SG.o

build/generalized_SG.a: build/generalized_SG.o build/distribution_table.o build/quad.o build/util.o
	ar rc build/generalized_SG.a build/generalized_SG.o build/quad.o build/util.o build/distribution_table.o

build/example: example.f90 build/generalized_SG.a
	${FC} ${FOPT} example.f90 -o build/example build/generalized_SG.a

build/C_interface.o: C_interface.f90 build/generalized_SG.o build/quad.o build/util.o build/distribution_table.o
	${FC} ${FOPT} -c C_interface.f90 -o build/C_interface.o

build/C_interface.a: build/C_interface.o build/generalized_SG.o build/quad.o build/util.o build/distribution_table.o
	ar rc build/C_interface.a build/C_interface.o build/generalized_SG.o build/quad.o build/util.o build/distribution_table.o

build/C_example: C_example.c ITHE.h build/C_interface.a
	${CC} ${COPT} C_example.c -o build/C_example build/C_interface.a ${CLIB}

clean:
	rm -rf build/*
