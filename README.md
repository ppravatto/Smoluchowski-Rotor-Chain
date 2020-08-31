# Smoluchowski-Rotor-Chain
Smoluchowski-Rotor-Chain is a simple C++ program devoted to the diffusive simulation of a chain of hydrodynamically independent 2D rotors interacting with a nearest-neighbour potential. The software solves the Symmetrized Fokker-Planck-Smoluchowski eigenvalue problem by expansion over a basis set.

## Requirements
* [gcc (g++)](https://gcc.gnu.org/) 9.3
* [CBLAS](http://www.netlib.org/blas/) 3.8
* [LAPACK](http://www.netlib.org/lapack/) 3.8
* [gnu-gsl](https://www.gnu.org/software/gsl/) 2.6
* [OpenMP](https://www.openmp.org/) (only for parallel version)

## Compiling
 A `Makefile` is provided to compile the code. The serial compilation can be linked against the system installation of the CBLAS library:
 ```
 make serial
 ```
 or against the CBLAS version provided by GNU-GSL:
 ```
 make serial-gsl
 ```
 A parallel (multi-threaded) version of the software can also be compiled linking against the OpenMP library either with the system CBLAS library:
 ```
 make parallel
 ```
 or against the one provided by GNU-GSL:
 ```
 make parallel-gsl
 ```
 The `make clean` command can be used to remove the executable file.
