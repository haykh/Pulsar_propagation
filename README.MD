### Libraries description
`lib` folder contains all the necessary libraries.
* `Nrutil.h` and `diffeqsolver.h` are from 'Numerical Recipes' textbook and contain everything to solve nonlinear ODE wirh rk4 routine.
* `integrator.h` is a custom one dimensional function integrator. `functions.h` contains custom functions and methods for vector analysis.
* All the physical functions involved are in `process_functions.h` and the right-hand side of the ODE is writen in `RHS.h`. `auxilary.h` is to store some of the auxilary routines used, such as finding the radiation point.
* `constants.h` is to store all the necessary constants.

### Running
* The input parameters are passed via the `input.parameters` file, from the same directory as `*.o` file.
* The compilation is done with `make all` in the root directory. The compiled file is being generated in the `/bin` directory.

_Note that `src/` is a temporary directory, created after `Makefile` is called. To adjust `main.cpp` and `pulsar_parameters.*` see `cases/*/src/`._

`.gitignore` file contains all the extensions and folders that we don't want to share to github, such as executables, images, outputs, python files, etc.

### Things to do
- Write density routine for any given ___B___-field.
- Make initial radiation point for any given magnetic field.
