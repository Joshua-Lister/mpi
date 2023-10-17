#pragma once
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <vector>
#include <mpi.h>

using namespace std;

/**
 The algorithm class which has all the neccessary methods to calculate 
 the 2D wave equation up to a user defined value in 'input.txt'.
 */
class algorithm
{
public:
	// Constructor
	algorithm(int imax, int jmax, string boundary_condition, int p, int id,
			  int x_max, int y_max, double tmax,
			  double dt_out, double t_out, int c, double dt);
	// Destructor
	~algorithm();

	//Organises given number of processors into a grid as square as possible
	void setup_domain(int &p, int &id);
	//Creates dynamic memory for 1D and 2D array and
	//maps a 1D array to a 2D array.
	void setup_continuous_array(double **&array_2d, double *&array_1d, int row_dim, int col_dim);
	// Frees dynamic memory of both 1D and 2D array
	void free_continuous_array(double **&array_2d, double *&array_1d);
	// Creates an initial splash condition
	void initial_condition();
	// Creates the relevant datatypes to recieve into ghost cells and send from actual cells
	void create_datatypes(double **&data, double *data_1D);
	//Fills send_data with mpi datatypes: _left, _top, _right and bottom
	//Fills recv_data according to boundary condition, if periodic then
	//edges are opposite.
	void direction();
	// Determines id of neighbouring processors depending on the boundary condition
	void calc_neighbours(string boundary_condition);
	// Sets up sends and recieves to neighbouring processors
	void do_comms();
	// Iterates from i = 2 to n - 2 while communications are being set up
	void do_interior_iteration();
	// Iterates over cells which include ghost cells in the stencil after communications are set up
	void do_iteration_ext_periodic_non_periodic_no_edge();
	// Iterates over all cells except ghost cells (this function is only used for comparison of
	// performance)
	void do_iteration_rubbish();
	//This assigns dirchlet boundary to the cells on the physical boundary which are
	// no longer zero due to receving columns & rows from other processors.
	void apply_dirchlet_bc_loop();
	// No physical boundary so if statements required. Faster to use boundary specific
	// function.
	void ext_iteration_periodic();
	// If number of rows = 1, only columns are sent with mpi, bottom and top row no mpi
	// comunication needed.
	void set_periodic_without_comms();
	// Applys dirchlet boundary condition (0 value on perimeter of new_grid, not the ghost cells!)
	void apply_dirchlet_bc();
	// Applys neumann boundary condition ∇u dot n = 0 (0 gradient in u direction normal to boundary)
	// Not set in the ghost cells !
	void apply_neumann_bc();
	// Runs the simulation with periodic boundaries if number of rows = 1
	void run_periodic_prime_processor();
	// Writes out local processor grid to a dat file at specific timesteps dt_out
	void grid_to_file(int out);
	// Runs the simulation with periodic boundaries
	void run_periodic();
	// Runs the simulation with dirchlet boundaries
	void run_dirchlet();
	// Runs the simulation with neumann boundaries
	void run_neumann();

	//! Intialises global grid rows and column dimensions to dummy value
	int imax = -1, jmax = -1;
	//! Initialises boundary condition to dummy value
	string boundary_condition = "";
	//! local row and column dimension after domain decomposition
	int i_loc_max = -1, j_loc_max = -1;
	//! local row and column dimension with ghost gells after domain decomposition
	int i_loc_w_ghost_max = -1, j_loc_w_ghost_max = -1;
	//! row and column index for division of processors grid
	int i_dom = -1, j_dom = -1;
	//! row and column dimension of division of processors grid
	int n = -1, m = -1;
	//! neighbours
	int _left = -1;
	int _top = -1;
	int _right = -1;
	int _bottom = -1;
	//! Contains id of neighbouring processors
	int neighbours[4];
	//! Contains mpi datatype for sending data
	MPI_Datatype send_data[4];
	//! Contains mpi datatype for receiving data into ghost cells
	MPI_Datatype recv_data[4];
	//! To track number of communications that have been setup in do_comms()
	int cnt = 0;
	//! Current grid 1D and 2D form
	double **grid = nullptr;
	double *grid1 = nullptr;
	//! New grid 1D and 2D form
	double **new_grid = nullptr;
	double *new_grid1 = nullptr;
	//! Old grid 1D and 2D form
	double **old_grid = nullptr;
	double *old_grid1 = nullptr;

	//! Physical parameters required to make initial splash condition
	double r_splash = 1.0;
	double x_splash = 3.0;
	double y_splash = 3.0;
	//! Increment from 0 to x_max and y_max by dx and dy respectively
	double dx = -1, dy = -1;
	//! Timestep initialised with dummy value
	double dt = -1;
	//! Constant c initialised with dummy value
	int c = -1;
	//! Starting time
	double t = 0;
	//! Iteration count
	int it = 0;
	//! Set to 1000 to avoid sorting problem in processing.py
	int out_cnt = 1000;
	//! Limit of range of x and y
	int x_max, y_max;
	//! break while loop when t_max is reached.
	double t_max;
	//! Time when grid file outputted.
	double t_out;
	//! The increment until next grid is outputted to file.
	double dt_out;
	// Tracks nonblocking communications
	MPI_Request *requests;
	//! MPI Datatypes for edges of grid both 1 cell in i,j and ghost cell edges.
	MPI_Datatype top_type, bottom_type, left_type, right_type,
		top_ghost_type, bottom_ghost_type, left_ghost_type, right_ghost_type;
	//! id of processor the program is being run on
	int id = -1;
	//! Number of total processors
	int p = -1;
};
