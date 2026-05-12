
# generalized_SG

Implementation of the Generalized Scharfetter-Gummel (SG) scheme using adaptive Tanh-Sinh quadrature. This project provides a flexible Fortran and C library for accurate and efficient computation of edge currents in drift-diffusion and related models, supporting arbitrary density of states and distribution functions.

## Features
- Generalized SG scheme for arbitrary density of states and distribution densities (Maxwell-Boltzmann/Fermi-Dirac)
- Adaptive Tanh-Sinh quadrature for robust integration
- Create lookup table for fast repeated evaluation of cumulative distribution function
- Fortran and C interface/example
- Internal parallelization with OpenMP

## Directory Structure
- `src/` - Fortran and C source files
	- `generalized_SG.f90` - Main Fortran module implementing the generalized SG scheme
	- `distribution.f90` - Example distribution functions
	- `fermi_dirac.f90`, `fermi_dirac_parameters.f90` - Fermi-Dirac integral support
	- `cheb_table.f90` - Chebyshev polynomial lookup table (requires Intel MKL)
	- `integral_table.f90` - Definite integral lookup table
	- `quad.F90` - Tanh-Sinh numerical quadrature routines (preprocessed with fypp)
	- `dual.F90` - Dual numbers for automatic differentiation (preprocessed with fypp)
	- `util.f90` - Utility functions
	- `mkl_dfti.f90` - Intel MKL Fast Fourier Transform interface
	- `ITHE.f90` - Fortran API for the C/Fortran interface library
	- `example.f90` - Fortran example program
	- `example.c` - C example program
- `include/` - Public C headers (`ITHE.h`, `distribution.h`)
- `build/` - Build artifacts (`*.o`, `*.mod`, static libraries, example binaries)
- `makefile` - Build instructions

## Build Instructions
Requirements:
- `gfortran` (GCC Fortran compiler)
- `gcc` (GCC C compiler)
- OpenMP support
- Intel MKL (for Chebyshev table acceleration via FFT)
- `fypp` (Fortran preprocessor)

To build the library and example programs, run:

```sh
make
```

This will produce the following in the `build/` directory:
- `generalized_SG.a` - Fortran static library with core SG scheme
- `ITHE.a` - Static library with C/Fortran interface
- `example` - Fortran example executable
- `example_c` - C example executable

The default compilation flags enable:
- Optimization level 2 (`-O2`) for performance
- OpenMP parallelization (`-fopenmp`)
- Execution stack support (`-z execstack`)
- Use double precision instead of single precision by default (`-fdefault-real-8`)

For debug builds with runtime checking and detailed diagnostics, uncomment the debug flags in the makefile (set `FOPT` to the commented line with `-O0 -g3 -ggdb -fcheck=all` etc.).

MKL linking can be switched between dynamic and static libraries in the makefile:
- Default (`FLIB`) uses dynamic MKL libraries (`-lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core`)
- An alternative static link line is provided in comments (using `-Wl,--start-group ...libmkl_*.a ... -Wl,--end-group`)

If you need fully static MKL linkage, replace the active `FLIB` line with the commented static variant and rebuild.

## Usage

### Fortran
You can link against `build/generalized_SG.a` and use the following modules in your code:
- `generalized_sg_m` - Main generalized SG scheme module
- `distribution_m` - Example distribution functions
- `fermi_dirac_m` - Fermi-Dirac integral
- `cheb_table_m` - Chebyshev table
- `integral_table_m` - Lookup table for definite integrals
- `quad_m` - Tanh-Sinh numerical quadrature routines
- `dual_m` -  Dual numbers for automatic differentiation
- `util_m` - Utility module

See [src/example.f90](src/example.f90) for usage details.

Run the Fortran example:
```sh
./build/example
```

### C
Link with `build/ITHE.a` and include `include/ITHE.h` and possibly `include/distribution.h` in your code. See [src/example.c](src/example.c) for usage details.

Run the C example:
```sh
./build/example_c
```

## Environment Setup

The build system requires the `MKLROOT` environment variable to be set to the Intel MKL installation directory. For example:

```sh
export MKLROOT=/opt/intel/mkl
```

## License
This project is licensed under the MIT License. See the `LICENSE` file for details. Example files may be under different licenses (see file headers).
