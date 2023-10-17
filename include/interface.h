#pragma once
#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
#include <mpi.h>
class interface
{
public:
	// constructor
	interface();
	// destructor
	~interface();
	//! Intialises global grid rows and column dimensions to dummy value
	int imax = -1;
	int jmax = -1;
	//! Increment from 0 to x_max and y_max by dx and dy respectively
	int x_max = -1;
	int y_max = -1;
	//! break when t_max reached
	double tmax = -1;
	//! increment until next grid is outputted to file
	double dt_out = -1;
	//! when grid file outputted
	double t_out = -1;
	//! physical parmeter
	int c = -1;
        //! timestep
	double dt = -1;
	//! boundary condition (Periodic, Neumann, Dirchlet)
	std::string boundary_condition = "";

private:
	// As private thus will not description will not be
	// included in doxygen index.html but can find in
	// interface.cpp
	void gatherdetails();
};
