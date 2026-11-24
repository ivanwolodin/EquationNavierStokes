# 2D Flow Simulation

A Fortran program for numerical simulation of a two-dimensional flow on a uniform 20 × 20 grid.

The program calculates the stream function and vorticity fields using finite difference methods.

## Features

- 2D uniform computational grid
- Finite difference scheme
- Iterative solution of the stream function equation
- Boundary condition implementation
- Export of results to `1.txt`

## Running

Compile the source file:

```bash
gfortran main.f90 -o simulation
```

Run the program:

```bash
./simulation
```

The results are saved to `1.txt`.