/*!
    *\author Joshua Piers Lister
    *\copyright MIT License
    *\mainpage Wave Equation 2D in parallel
    *\section intro_sec Introduction
    * This is an implementation of the two dimensional(2D) wave equation in c++ using parallel methods 
      with Message Passing Interface(MPI). The program contains the algorithm class in /src which contains 
      all of the methods required to calculate the wave equation and a main file which implements a boundary 
      condition based on the user input to the text file 'input.txt'. It also contains the interface class 
      which reads in 'input.txt' assigns values to reaction condition variables. The file processing.py reads 
      in the data outputted by the function grid_to_file in algorithm class and plots a contourf animation. 
      It is a standalone program and has no dependicies.
      More details can be found under Related Pages tab and click acse-6-mpi-coursework-acse-jpl20.
    *\subsection compile_sec Compilation
    * Enter into the bash command line 'make', after successful build enter mpiexec -n ? run_parallel where ?
      is the number of processors you wish to use e.g. mpiexec -n 2 run_parallel.
    */
#include "algo.h"
#include "interface.h"
#include <chrono>
#include <time.h>
#include <iomanip>

#define DO_TIMING

using namespace std;

int id, p, tag_num = 1;

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
#ifdef DO_TIMING
    MPI_Barrier(MPI_COMM_WORLD);
    auto start1 = chrono::high_resolution_clock::now();
#endif

    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    interface collect;


    auto *run = new algorithm(collect.imax, collect.jmax, collect.boundary_condition, p, id, collect.x_max, collect.y_max,
                              collect.tmax, collect.dt_out, collect.t_out, collect.c, collect.dt);

    if (collect.boundary_condition == "Periodic" && run->m == 1)
        run->run_periodic_prime_processor();
    else if (collect.boundary_condition == "Dirchlet")
        run->run_dirchlet();
    else if (collect.boundary_condition == "Neumann")
        run->run_neumann();
    else if (collect.boundary_condition == "Periodic")
        run->run_periodic();
    else
    {
        MPI_Finalize();
        cerr << "No valid boundary condition entered\n";
        exit(0);
    }
#ifdef DO_TIMING
    MPI_Barrier(MPI_COMM_WORLD);
    auto end1 = chrono::high_resolution_clock::now();
    if (id == 0)
    {
        std::chrono::duration<double> elapsed = end1 - start1;
        cout << setprecision(5);
        cout << "The code took " << elapsed.count() << "s to run on " << p
             << " processors with boundary condition " << collect.boundary_condition << endl;
        std::stringstream ffname;
        std::fstream f2;
        ffname << "./timings/" << collect.boundary_condition << "_good_algorithm.dat";
        f2.open(ffname.str().c_str(), fstream::app);
        f2 << p << "\t" << elapsed.count() << "\n";
    }
#endif
    MPI_Finalize();
return 0;
}
