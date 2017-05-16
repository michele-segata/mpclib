# MPCLib

`MPCLib` is a C++ library that solves Model Predictive Control problems. It can
also be used to solve strictly convex quadratic minimization problems by using
the embedded [QuadProg++](https://github.com/liuq/QuadProgpp) library.

The library defines a very simple interface that can be used to enable control
and output constraints, slack variables, minimization weights, etc. In turn, the
MPC interface uses the `QuadraticProblem` and the `Constraint` classes to easily
define constraints and add them to the minimization problem.

Please be aware that this software is still in its early life. It works fine but
it is not complete. The software still doesn't provide useful errors or
exceptions when something is misconfigured. The project will be updated over
time. Feel free to contribute.

## Building and running the sample application

The project is `cmake` based. To build it simply type

```
mkdir build
cd build
cmake ..
make
```

To run some tests and generate some plots simply run `test-mpc.sh` in the main
folder. The script runs 6 different tests producing `csv` files as output. If
you have [R](https://www.r-project.org/) installed, the script will also plot
the results. For more formal information on the problem the library solves and
on the examples please refer to the PDF in the `doc` folder.

## Included software

The library includes other open source software (please see the LICENSE file):

 * [QuadProg++](https://github.com/liuq/QuadProgpp): this is the library that
   solves the quadratic optimization problem;
 * [TCLAP](https://github.com/eile/tclap): this is a command line arguments
   parsing library that is used in the sample program;

## Contribution

Copyright (C) 2017 Michele Segata
