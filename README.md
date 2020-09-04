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

## Running calculations
In order to perform a calulation an input file must be provided as inline argument to the `main.exe` executable.
```
./main.exe input.txt
```
The program can be executed in background and the output redirected to an `output.txt` file using the syntax:
```
./main.exe input.txt > output.txt &
```
If a parallel version of the program is compiled the number of threads can be selected through the OpenMP environment variable `OMP_NUM_THREADS`.

### Input file syntax
The command syntax for the program must be always structurated as follows:
```
# KEYWORD
parameters
```
The line marked by `#` identify a command line. The line should contain the marker `#` followed by a single white space and a keyword. The full list of keyword is available in the next section. The first line after the keyword must contain the parameters either as a single number (if only one is required) or as a list of numbers separated by a `,`. Floating point numbers can be expressed in exponential form (e.g. `6.02e-23`). The `.` must be used as decimal point.

Every input file should contain the number of rotors to be considered. This can be specified as:
```
# NROT
number of rotors
```
Please note how `# NROT` represent the number of rotors and not the number of dihedral angles.

### Input file keywords
Keyword | Description | Number of parameters | Type | Default value
 :--- | :---: | :---: | :---: | :---:
`NROT` | Number of rotors in the chain | 1 | `integer` | -
`NUM-MIN` | Number of potential minimum | `NROT-1` | `integer list` | 1
`ISOLATED-BASIS` | Number of single dihedral basis functions | `NROT-1` | `integer list` | 100
`COUPLED-BASIS` | Number of single dihedral eigenfunctions to be used in the multiple rotor expansion | `NROT-1` | `integer list` | 10
`BARRIER` | Potential energy barrier for each dihedral | `NROT-1` | `double list` | 1.
`DIFFUSION-COEFF` | Diffusion coefficient for each rotor | `NROT` | `double list` | 1.
`NPT-INTEG-SINGLE` | Maximum number of integration points for the single dihedral | 1 | `int` | 10000
`ABS-INTEG-SINGLE` | Absolute integrator error threshold for the single dihedral | 1 | `double` | 1e-10
`REL-INTEG-SINGLE` | Relative integrator error threshold for the single dihedral | 1 | `double` | 1e-10
`QAG-KEY-SINGLE` | GSL key for selecting the QAG integrator order for the single dihedral | 1 | `int` | 6
`NPT-INTEG-COUPLED` | Maximum number of integration points for the rotor chain | 1 | `int` | 10000
`ABS-INTEG-COUPLED` | Absolute integrator error threshold for the rotor chain | 1 | `double` | 1e-10
`REL-INTEG-COUPLED` | Relative integrator error threshold for the rotor chain | 1 | `double` | 1e-10
`QAG-KEY-COUPLED` | GSL key for selecting the QAG integrator order for the rotor chain | 1 | `int` | 6
