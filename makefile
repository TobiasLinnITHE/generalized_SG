# SPDX-License-Identifier: MIT
# This file is licensed under the MIT License.
# See https://opensource.org/licenses/MIT

FYPP = fypp
FC   = gfortran
# FOPT = -ffree-line-length-none -fopenmp -fdefault-real-8 -z execstack -J build/ -O0 -g3 -ggdb -fcheck=all -fbacktrace -finit-integer=-9999 -finit-real=-inf -ffpe-trap=invalid,zero,overflow
FOPT = -ffree-line-length-none -fopenmp -fdefault-real-8 -z execstack -J build/ -O2

MKLOPT     = -I"${MKLROOT}/include"
# FLIB = -L${MKLROOT}/lib -Wl,--start-group ${MKLROOT}/lib/libmkl_gf_lp64.a ${MKLROOT}/lib/libmkl_gnu_thread.a ${MKLROOT}/lib/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl
FLIB = -L${MKLROOT}/lib -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl

CC   = gcc
COPT = -z execstack -O2 -I"include"
CLIB = -lgfortran -lm -lquadmath ${FLIB}

AR    = ar
AROPT = rc



.PHONY: all clean



all: build build/generalized_SG.a build/ITHE.a build/example build/example_c




build:
	mkdir -p build

build/util.o: src/util.f90
	${FC} ${FOPT} -c src/util.f90 -o build/util.o

src/dual.f90: src/dual.F90
	${FYPP} src/dual.F90 > src/dual.f90

build/dual.o: src/dual.f90 build/util.o
	${FC} ${FOPT} -c src/dual.f90 -o build/dual.o

src/quad.f90: src/quad.F90
	${FYPP} src/quad.F90 > src/quad.f90

build/quad.o: src/quad.f90
	${FC} ${FOPT} -c src/quad.f90 -o build/quad.o

build/integral_table.o: src/integral_table.f90 build/util.o
	${FC} ${FOPT} -c src/integral_table.f90 -o build/integral_table.o

build/cheb_table.o: src/cheb_table.f90 build/util.o
	${FC} ${FOPT} ${MKLOPT} -c src/cheb_table.f90 -o build/cheb_table.o

build/generalized_SG.o: src/generalized_SG.f90 build/quad.o build/dual.o build/util.o
	${FC} ${FOPT} -c src/generalized_SG.f90 -o build/generalized_SG.o

build/generalized_SG.a: build/generalized_SG.o build/cheb_table.o build/integral_table.o build/quad.o build/dual.o build/util.o
	${AR} ${AROPT} build/generalized_SG.a build/generalized_SG.o build/cheb_table.o build/integral_table.o build/quad.o build/dual.o build/util.o



build/fermi_dirac_parameters.o: src/fermi_dirac_parameters.f90 build/util.o
	${FC} ${FOPT} -c src/fermi_dirac_parameters.f90 -o build/fermi_dirac_parameters.o

build/fermi_dirac.o: src/fermi_dirac.f90 build/fermi_dirac_parameters.o build/util.o
	${FC} ${FOPT} -c src/fermi_dirac.f90 -o build/fermi_dirac.o

build/distribution.o: src/distribution.f90 build/fermi_dirac.o build/integral_table.o build/cheb_table.o
	${FC} ${FOPT} -c src/distribution.f90 -o build/distribution.o

build/example: src/example.f90 build/distribution.o build/fermi_dirac.o build/fermi_dirac_parameters.o build/generalized_SG.a
	${FC} ${FOPT} ${MKLOPT} src/example.f90 -o build/example build/distribution.o build/fermi_dirac.o build/fermi_dirac_parameters.o build/generalized_SG.a ${FLIB}



build/ITHE.o: src/ITHE.f90 build/generalized_SG.a
	${FC} ${FOPT} -c src/ITHE.f90 -o build/ITHE.o

build/ITHE.a: build/ITHE.o build/generalized_SG.o build/cheb_table.o build/integral_table.o build/quad.o build/dual.o build/util.o
	${AR} ${AROPT} build/ITHE.a build/ITHE.o build/generalized_SG.o build/cheb_table.o build/integral_table.o build/quad.o build/dual.o build/util.o



build/example_c: src/example.c include/ITHE.h include/distribution.h build/distribution.o build/fermi_dirac.o build/fermi_dirac_parameters.o build/ITHE.a
	${CC} ${COPT} src/example.c -o build/example_c build/distribution.o build/fermi_dirac.o build/fermi_dirac_parameters.o build/ITHE.a ${CLIB}



clean:
	rm -rf build/*
	rm src/dual.f90
	rm src/quad.f90
