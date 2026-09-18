# Exact Diagonalization Code for BHZ model

### Features

- Flexible with respect to boundary conditions (PBC or OBC) along both directions x and y
- Works in both canonical and grand-canonical ensemble
- Computes local charge and magnetic properties
- Calculates topological band structures, spin currents, edge state wave functions
 

### Requirements

- CMake
- g++
- C++11 or newer
- LAPACK and BLAS

### Details and Compilation

To compile the code do:
```bash
make -f Makefile
```
or
```bash
make
```

To run the code do:
```bash
./ED_exe input.inp
```
