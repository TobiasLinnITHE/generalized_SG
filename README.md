
# generalized_SG

Implementation of the Generalized Scharfetter-Gummel (SG) scheme using adaptive Tanh-Sinh quadrature. This project provides a flexible Fortran and C library for accurate and efficient computation of edge currents in drift-diffusion and related models, supporting arbitrary density of states and distribution functions.

## Features
- Generalized SG scheme for arbitrary density of states and distribution densities (Maxwell-Boltzmann/Fermi-Dirac)
- Adaptive Tanh-Sinh quadrature for robust integration
- Create lookup table for fast repeated evaluation of cumulative distribution function
- Fortran and C interface/example
- Internal parallelization with OpenMP

## Directory Structure
- `generalized_SG.f90` — Main Fortran module implementing the generalized SG scheme
- `distribution_table.f90` — Lookup table for cumulative distribution functions
- `quad.f90` — Tanh-Sinh quadrature routines
- `util.f90` — Utility functions
- `C_interface.f90`, `C_example.c`, `ITHE.h` — C interface and example usage
- `example.f90` — Fortran example program
- `makefile` — Build instructions

## Build Instructions
Requirements: `gfortran`, `gcc`, OpenMP support

To build the library and example programs, run:

```sh
make -j4
```

This will produce the following in the `build/` directory:
- `generalized_SG.a` — Fortran static library
- `C_interface.a` — C interface static library
- `example` — Fortran example executable
- `C_example` — C example executable

## Usage

### Fortran
You can link against `build/generalized_SG.a` and use the generalized_sg_m (`build/generalized_sg_m.mod`) and distribution_table_m (`build/distribution_table_m.mod`) modules in your code. See `example.f90` for more details.

Run the Fortran example:
```sh
./build/example
```

### C
Link with `build/C_interface.a` (C) in your code, use the `ITHE.h` header file. See `C_example.c` for more details.

Run the C example:
```sh
./build/C_example
```

## License
This project is licensed under the MIT License. See the `LICENSE` file for details. Example files may be under different licenses (see file headers).
