# Parallelised implementation of 2D wave equation using MPI 

This is an implementation of the two dimensional(2D) wave equation in c++ using parallel methods with Message Passing Interface(MPI). The program contains the algorithm class in /src which contains all of the methods required to calculate the wave equation and a main file which implements a boundary condition based on the user input to the text file 'input.txt'. It also contains the interface class which reads in 'input.txt' assigns values to reaction condition variables. The file processing.py reads in the data outputted by the function grid_to_file in algorithm class ; plots ands saves an animation. It is a standalone program and has no dependicies. Visit docs/html/index.html for information on code documentation.

## Requirements

This library requires a compiler that supports C++17 language standard and has the library <mpi.h> installed. To run processing.py to produce the animation the user requires either a conda environment or the following modules for python installed: pandas, matplotlib and glob. 

## Algo

Algo contains the algorithm class which has all the neccessary methods to calculate the 2D wave equation up to a user defined value in 'input.txt'.

### Notable Properties

- `i_loc_max`(`int`): grid rows dimension.
- `j_loc_max`(`int`) :  grid columns dimension.
- `i_loc_w_ghost_max`(`int`): grid rows dimension plus ghost/halo cells.
- `j_loc_w_ghost_max`(`int`) : grid dimension columns plus ghost/halo cells.
- `_left`(`int`): left neighbour id.
- `_top`(`int`):  top neighbour id.
- `_right`(`int`):  right neighbour id.
- `_bottom`(`int`): bottom neighbour id.
- `r_splash`(`double`): initial condition r splash value.
- `x_splash`(`double`): initial condition x splash value.
- `y_splash`(`double`): initial condition y splash value.

### Notable Methods

- `void setup_domain(int &p, int &id)`: Splits up a number of processors p into an optimal configuration.
- `void setup_continuous_array(double**& array_2d, double*& array_1d, int r, int c)`: Maps a 1D array to a 2D one.
- `void initial_condition()`: Creates the initial splash condition based on the global i and j values.
- `void create_datatypes(double** &data, double* data_1D)`: Creates datatypes for all four directions for both ghost cells and cells one in from the ghost cells.
- `void calc_neighbours(string boundary_condition)`:  Calculates the neighbouring ids depending on the boundary_condition value.
- `void do_comms()`:  Sets up sends and receives.
- `do_interior_iteration()`:  Calculates the interior nodes which don't depend on any ghost cells whilst do_comms() is still in progress.
- `void do_iteration_ext_periodic_non_periodic_no_edge`: Calculates exterior nodes which do depend on ghost cells. Is only called after do_comms() is complete.

## Interface

Interface contains a method to read in a text file and assign the values read in to the appropriate variables.

## Running the program

On Linux based OS the user can compile the program in the terminal by entering `make` and then `mpiexec -n ? run_parallel` where `?` is to be replaced by the number of nodes you wish to run the program on.

## Installation

- Conda installation for linux - https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html
- Conda installation for windows - https://conda.io/projects/conda/en/latest/user-guide/install/windows.html
- Conda installation for mac - https://docs.conda.io/projects/continuumio-conda/en/latest/user-guide/install/macos.html
###

To have the neccessary modules required to run processing.py without conda  enter `pip install -r requirements.txt`.

## Documentation

Extensive code documentation can be found in docs/html/index.html. You will need to open it in your local web browser e.g. google-chrome docs/html/index.html

## Licence 

[MIT](https://choosealicense.com/licenses/mit/)

