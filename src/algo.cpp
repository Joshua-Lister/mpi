#define _USE_MATH_DEFINES
#include "algo.h"
#include <chrono>

//!Constructor
/**
    Assigns the variables required to run the simulation
    in parallel for dirchlet, neumann and periodic boundary
    conditions and to write grid to file at specified times.
    Automatically calls setup_domain to arrange processors into
    optimal grid.
*/
algorithm::algorithm(int imax, int jmax, string boundary_condition, int p, int id, int x_max,
                     int y_max, double t_max, double dt_out, double t_out, int c, double dt) : imax(imax),
                                                                                               jmax(jmax), boundary_condition(boundary_condition), p(p), id(id), x_max(x_max),
                                                                                               y_max(y_max), t_max(t_max), dt_out(dt_out), t_out(t_out), c(c), dt(dt)
{
    setup_domain(p, id);
}

//! Destructor
/*!
    The destructor is called at the end of the program
    and free alls mpi types and the dynamic memory used
    for the 1D and 2D arrays for grid, new_grid and old_grid.
*/
algorithm::~algorithm()
{
    MPI_Type_free(&bottom_type);
    MPI_Type_free(&left_type);
    MPI_Type_free(&right_type);
    MPI_Type_free(&top_type);
    MPI_Type_free(&top_ghost_type);
    MPI_Type_free(&bottom_ghost_type);
    MPI_Type_free(&left_ghost_type);
    MPI_Type_free(&right_ghost_type);

    free_continuous_array(grid, grid1);
    free_continuous_array(new_grid, new_grid1);
    free_continuous_array(old_grid, old_grid1);
}

//! setup_domain
/*!
    Organises processors into optimal grid layout.
*/
void algorithm::setup_domain(int &p, int &id)
{
    // m = rows
    // n = cols
    m = 1;
    n = p;

    for (int i = 2; i < p; i++)
    {
        if (p % i == 0)
        {
            if (abs(i - p / i) < abs(m - n))
            {
                m = i;
                n = p / i;
            }
        }
    }

    i_dom = id / n;
    j_dom = id % n;

    i_loc_max = this->imax / m;
    j_loc_max = this->jmax / n;
    i_loc_w_ghost_max = (this->i_loc_max + 2);
    j_loc_w_ghost_max = (this->j_loc_max + 2);
}

//! initial_condition
/*!
    Initial splash condition calculated and assigned to 
    relevant grid entries.
*/
void algorithm::initial_condition()
{
    dx = x_max / ((double)imax - 1);
    dy = y_max / ((double)jmax - 1);
    for (int i = 2; i < i_loc_max - 1; i++)
        for (int j = 2; j < i_loc_max - 1; j++)
        {
            double x = dx * (i - 1) + (i_dom * i_loc_max);
            double y = dy * (j - 1) + (j_dom * j_loc_max);
            double dist = sqrt((x - x_splash) * (x - x_splash) + (y - y_splash) * (y - y_splash));
            if (dist < r_splash)
            {
                double h = 5.0 * (cos(dist / r_splash * M_PI) + 1.0);

                grid[i][j] = h;
                old_grid[i][j] = h;
            }
        }
}

//! calc_neighbours
/*!
    Calculates the id of the neighbouring processors
    depending on the boundary condition.
*/
void algorithm::calc_neighbours(string boundary_condition)
{
    
    _left = id - 1;
    _right = id + 1;
    _top = id - n;
    _bottom = id + n;
    if (boundary_condition == "Periodic" && m != 1)
    {
        if (j_dom == 0)
            _left = id + (n - 1);

        if (j_dom == n - 1)
            _right = id - (n - 1);

        if (i_dom == 0)
            _top = id + ((m - 1) * n);

        if (i_dom == m - 1)
            _bottom = id - ((m - 1) * n);
    }
    else if (boundary_condition == "Periodic" && m == 1)
    {
        if (j_dom == 0)
            _left = id + (n - 1);

        if (j_dom == n - 1)
            _right = id - (n - 1);

        _top = MPI_PROC_NULL;
        _bottom = MPI_PROC_NULL;
    }
    else
    {
        if (j_dom == 0)
            _left = MPI_PROC_NULL;

        if (j_dom == n - 1)
            _right = MPI_PROC_NULL;

        if (i_dom == 0)
            _top = MPI_PROC_NULL;

        if (i_dom == m - 1)
            _bottom = MPI_PROC_NULL;
    }

    neighbours[0] = _left;
    neighbours[1] = _top;
    neighbours[2] = _right;
    neighbours[3] = _bottom;
}

//! setup_continous_array
/*!
    Allocates dynamic memory (heap) to entry parameter
    array_1d of size i_loc_w_ghost_max * j_loc_w_ghost_max.
    It also allocates dynamic memory to entry parameter array_2d
    but of size j_loc_w_ghost and contains the memory address of
    iterable i * col_dim.
*/
void algorithm::setup_continuous_array(double **&array_2d, double *&array_1d, int row_dim, int col_dim)
{
    array_1d = new double[row_dim * col_dim];
    array_2d = new double *[row_dim];

    for (int i = 0; i < row_dim; i++)
        array_2d[i] = &array_1d[i * col_dim];
}

//! free_continous_array
/*!
    Frees dynamic memory (heap) associated with entry
    parameters array_2d and array_1d.
*/
void algorithm::free_continuous_array(double **&array_2d, double *&array_1d)
{
    delete[] array_1d;
    delete[] array_2d;
}

//! create_datatypes
/*!
    Creates datatypes to send from left, top, right and bottom and
    datatypes to recieve into left_ghost, top_ghost, right_ghost
    and bottom_ghost
*/
void algorithm::create_datatypes(double **&data, double *data_1D)
{
    // auto *block_lengths = new int[i_loc_max];
    // auto *typelist = new MPI_Datatype[i_loc_max];
    // auto *addresses = new MPI_Aint[i_loc_max];
    vector<int> block_lengths(i_loc_max);
    vector<MPI_Datatype> typelist(i_loc_max);
    vector<MPI_Aint> addresses(i_loc_max);

    MPI_Aint add_start;

    // RIGHT
    for (int i = 1; i <= i_loc_max; i++)
    {
        block_lengths[i - 1] = 1;
        typelist[i - 1] = MPI_DOUBLE;
        MPI_Aint temp_add;
        MPI_Get_address(&data[i][j_loc_max], &temp_add);
        addresses[i - 1] = temp_add;
    }
    MPI_Get_address(data_1D, &add_start);
    for (int i = 0; i < i_loc_max; i++)
        addresses[i] = addresses[i] - add_start;
    MPI_Type_create_struct(i_loc_max, block_lengths.data(), addresses.data(), typelist.data(), &right_type);
    MPI_Type_commit(&right_type);

    // LEFT
    for (int i = 1; i <= i_loc_max; i++)
    {
        MPI_Aint temp_add;
        MPI_Get_address(&data[i][1], &temp_add);
        addresses[i - 1] = temp_add;
    }
    for (int i = 0; i < i_loc_max; i++)
        addresses[i] = addresses[i] - add_start;
    MPI_Type_create_struct(i_loc_max, block_lengths.data(), addresses.data(), typelist.data(), &left_type);
    MPI_Type_commit(&left_type);

    //TOP
    int block_length = j_loc_max;
    MPI_Datatype typeval = MPI_DOUBLE;
    MPI_Aint address;
    MPI_Get_address(&data[1][1], &address);
    address = address - add_start;
    MPI_Type_create_struct(1, &block_length, &address, &typeval, &top_type);
    MPI_Type_commit(&top_type);

    //BOTTOM
    MPI_Get_address(&data[i_loc_max][1], &address);
    address = address - add_start;
    MPI_Type_create_struct(1, &block_length, &address, &typeval, &bottom_type);
    MPI_Type_commit(&bottom_type);

    //TOP GHOST
    MPI_Get_address(&data[0][1], &address);
    address = address - add_start;
    MPI_Type_create_struct(1, &block_length, &address, &typeval, &top_ghost_type);
    MPI_Type_commit(&top_ghost_type);

    //BOTTOM GHOST
    MPI_Get_address(&data[j_loc_w_ghost_max - 1][1], &address);
    address = address - add_start;
    MPI_Type_create_struct(1, &block_length, &address, &typeval, &bottom_ghost_type);
    MPI_Type_commit(&bottom_ghost_type);

    // RIGHT GHOST
    for (int i = 1; i <= i_loc_max; i++)
    {
        MPI_Aint temp_add;
        MPI_Get_address(&data[i][j_loc_w_ghost_max - 1], &temp_add);
        addresses[i - 1] = temp_add;
    }
    for (int i = 0; i < i_loc_max; i++)
        addresses[i] = addresses[i] - add_start;
    MPI_Type_create_struct(i_loc_max, block_lengths.data(), addresses.data(), typelist.data(), &right_ghost_type);
    MPI_Type_commit(&right_ghost_type);

    // LEFT GHOST
    for (int i = 1; i <= i_loc_max; i++)
    {
        MPI_Aint temp_add;
        MPI_Get_address(&data[i][0], &temp_add);
        addresses[i - 1] = temp_add;
    }
    for (int i = 0; i < i_loc_max; i++)
        addresses[i] = addresses[i] - add_start;
    MPI_Type_create_struct(i_loc_max, block_lengths.data(), addresses.data(), typelist.data(), &left_ghost_type);
    MPI_Type_commit(&left_ghost_type);

    // delete[] block_lengths;
    // delete[] typelist;
    // delete[] addresses;
}

//! direction
/*!
    Assigns the directions to send_data. If the boundary condition
    is periodic then the recv_data entries edges are flipped while if
    the boundary condition is neumann or dirchlet then recv_data entries
    remain the same direction as send_data.
*/
void algorithm::direction()
{
    send_data[0] = left_type;
    send_data[1] = top_type;
    send_data[2] = right_type;
    send_data[3] = bottom_type;

    if (boundary_condition == "Periodic")
    {
        recv_data[0] = right_ghost_type;
        recv_data[1] = bottom_ghost_type;
        recv_data[2] = left_ghost_type;
        recv_data[3] = top_ghost_type;
    }
    else
    {
        recv_data[0] = left_ghost_type;
        recv_data[1] = top_ghost_type;
        recv_data[2] = right_ghost_type;
        recv_data[3] = bottom_ghost_type;
    }
}

//! do_comms
/*!
    Loops through all the neighbouring ids of the processors
    and sets up the appropriate sends and receives. The variable
    cnt is is incremented such that MPI_Waitall can track the progress
    of the communications. If the boundary condition is not periodic,
    at edge cases, no sends or receives occur if the neighbour value is
    MPI_PROC_NULL.
*/
void algorithm::do_comms()
{
    requests = new MPI_Request[8];
    cnt = 0;
    for (int i = 0; i < 4; i++)
    {
        MPI_Irecv(grid1, 1, recv_data[i], neighbours[i], 0, MPI_COMM_WORLD, &requests[cnt]);
        cnt++;
    }
    for (int i = 0; i < 4; i++)
    {
        MPI_Isend(grid1, 1, send_data[i], neighbours[i], 0, MPI_COMM_WORLD, &requests[cnt]);
        cnt++;
    }
    // Do interior iteration while communications is being setup to reduce
    // computation time
    do_interior_iteration();
    // Wait until all sends & receives have finished until moving on.
    MPI_Waitall(cnt, requests, MPI_STATUS_IGNORE);
    delete[] requests;
}

//! do_interior_iteration
/*!
    Calculates interior nodes which do not require ghost cell values.
*/
void algorithm::do_interior_iteration()
{
    for (int i = 2; i < this->i_loc_w_ghost_max - 2; i++)
        for (int j = 2; j < this->j_loc_w_ghost_max - 2; j++)
            new_grid[i][j] = (dt * dt * c * c) * ((grid[i + 1][j] - 2.0 * grid[i][j] + grid[i - 1][j]) / (dx * dx) + (grid[i][j + 1] - 2.0 * grid[i][j] + grid[i][j - 1]) / (dy * dy)) + 2.0 * grid[i][j] - old_grid[i][j];
}

//! do_iteration_ext_periodic_non_periodic_no_edge
/*!
    Calculates exterior nodes which do require ghost cell values.
    This function is only called after all communications in do_comms()
    have finished successfully.
*/
void algorithm::do_iteration_ext_periodic_non_periodic_no_edge()
{

    // LEFT EDGE BUT ONE
    if (neighbours[0] != MPI_PROC_NULL)
    {
        int j = 1;
        for (int i = 1; i < this->i_loc_w_ghost_max - 1; i++)
            new_grid[i][j] = (dt * dt * c * c) * ((grid[i + 1][j] - 2.0 * grid[i][j] + grid[i - 1][j]) / (dx * dx) + (grid[i][j + 1] - 2.0 * grid[i][j] + grid[i][j - 1]) / (dy * dy)) + 2.0 * grid[i][j] - old_grid[i][j];
    }
    if (neighbours[2] != MPI_PROC_NULL)
    {
        // RIGHT EDGE BUT ONE
        int j = j_loc_w_ghost_max - 2;
        for (int i = 1; i < this->i_loc_w_ghost_max - 1; i++)
            new_grid[i][j] = (dt * dt * c * c) * ((grid[i + 1][j] - 2.0 * grid[i][j] + grid[i - 1][j]) / (dx * dx) + (grid[i][j + 1] - 2.0 * grid[i][j] + grid[i][j - 1]) / (dy * dy)) + 2.0 * grid[i][j] - old_grid[i][j];
    }
    if (neighbours[1] != MPI_PROC_NULL)
    {
        // TOP BUT ONE
        int i = 1;
        for (int j = 1; j < this->j_loc_w_ghost_max - 1; j++)
            new_grid[i][j] = (dt * dt * c * c) * ((grid[i + 1][j] - 2.0 * grid[i][j] + grid[i - 1][j]) / (dx * dx) + (grid[i][j + 1] - 2.0 * grid[i][j] + grid[i][j - 1]) / (dy * dy)) + 2.0 * grid[i][j] - old_grid[i][j];
    }
    if (neighbours[3] != MPI_PROC_NULL)
    {
        // BOTTTOM BUT ONE
        int i = i_loc_w_ghost_max - 2;
        for (int j = 1; j < this->j_loc_w_ghost_max - 1; j++)
            new_grid[i][j] = (dt * dt * c * c) * ((grid[i + 1][j] - 2.0 * grid[i][j] + grid[i - 1][j]) / (dx * dx) + (grid[i][j + 1] - 2.0 * grid[i][j] + grid[i][j - 1]) / (dy * dy)) + 2.0 * grid[i][j] - old_grid[i][j];
    }
}
void algorithm::ext_iteration_periodic()
{
    // LEFT EDGE BUT ONE
    int j = 1;
    for (int i = 1; i < this->i_loc_w_ghost_max - 1; i++)
        new_grid[i][j] = (dt * dt * c * c) * ((grid[i + 1][j] - 2.0 * grid[i][j] + grid[i - 1][j]) / (dx * dx) + (grid[i][j + 1] - 2.0 * grid[i][j] + grid[i][j - 1]) / (dy * dy)) + 2.0 * grid[i][j] - old_grid[i][j];

    // RIGHT EDGE BUT ONE
    j = j_loc_w_ghost_max - 2;
    for (int i = 1; i < this->i_loc_w_ghost_max - 1; i++)
        new_grid[i][j] = (dt * dt * c * c) * ((grid[i + 1][j] - 2.0 * grid[i][j] + grid[i - 1][j]) / (dx * dx) + (grid[i][j + 1] - 2.0 * grid[i][j] + grid[i][j - 1]) / (dy * dy)) + 2.0 * grid[i][j] - old_grid[i][j];

    // TOP BUT ONE
    int i = 1;
    for (int j = 1; j < this->j_loc_w_ghost_max - 1; j++)
        new_grid[i][j] = (dt * dt * c * c) * ((grid[i + 1][j] - 2.0 * grid[i][j] + grid[i - 1][j]) / (dx * dx) + (grid[i][j + 1] - 2.0 * grid[i][j] + grid[i][j - 1]) / (dy * dy)) + 2.0 * grid[i][j] - old_grid[i][j];

    // BOTTTOM BUT ONE
    i = i_loc_w_ghost_max - 2;
    for (int j = 1; j < this->j_loc_w_ghost_max - 1; j++)
        new_grid[i][j] = (dt * dt * c * c) * ((grid[i + 1][j] - 2.0 * grid[i][j] + grid[i - 1][j]) / (dx * dx) + (grid[i][j + 1] - 2.0 * grid[i][j] + grid[i][j - 1]) / (dy * dy)) + 2.0 * grid[i][j] - old_grid[i][j];
}

//! do_iteration_rubbish
/*!
    Calculates all nodes. This function is only used to compare
    performance against setting up comms, doing interior nodes then
    doing exterior nodes after comms have finished.
*/
void algorithm::do_iteration_rubbish()
{
    for (int i = 1; i <= i_loc_max; i++)
        for (int j = 1; j <= j_loc_max; j++)
            new_grid[i][j] = (dt * dt * c * c) * ((grid[i + 1][j] - 2.0 * grid[i][j] + grid[i - 1][j]) / (dx * dx) + (grid[i][j + 1] - 2.0 * grid[i][j] + grid[i][j - 1]) / (dy * dy)) + 2.0 * grid[i][j] - old_grid[i][j];
}

//! apply_dirchlet_bc
/*!
    Checks if neighbours are out of bounds,
    if so the physical boundary is set to 0.
*/
void algorithm::apply_dirchlet_bc()
{
    // LEFT
    if (neighbours[0] == MPI_PROC_NULL)
    {
        for (int i = 1; i <= i_loc_max; i++)
            new_grid[i][1] = 0;
    }
    // TOP
    if (neighbours[1] == MPI_PROC_NULL)
    {
        for (int j = 1; j <= j_loc_max; j++)
            new_grid[1][j] = 0;
    }
    // RIGHT
    if (neighbours[2] == MPI_PROC_NULL)
    {
        for (int i = 1; i <= i_loc_max; i++)
            new_grid[i][j_loc_w_ghost_max - 2] = 0;
    }
    // BOTTOM
    if (neighbours[3] == MPI_PROC_NULL)
    {
        for (int j = 1; j <= j_loc_max; j++)
            new_grid[i_loc_w_ghost_max - 2][j] = 0;
    }
}

//!apply_dirchlet_bc_loop
/*
 * Applys the dirchlet boundary condition 
 * only to the grid cells on the physical boundary
 * whose value changes from 0 due to exchange of
 * rows & columns from neighbouring processors.
 */
void algorithm::apply_dirchlet_bc_loop()
{
    // LEFT
    if (neighbours[0] == MPI_PROC_NULL)
    {
        new_grid[i_loc_max][1] = 0;
        new_grid[1][1] = 0;
    }
    // TOP
    if (neighbours[1] == MPI_PROC_NULL)
    {
        new_grid[1][j_loc_max] = 0;
        new_grid[1][1] = 0;
    }
    // RIGHT
    if (neighbours[2] == MPI_PROC_NULL)
    {
        new_grid[i_loc_max][j_loc_max] = 0;
        new_grid[1][j_loc_max] = 0;
    }
    // BOTTOM
    if (neighbours[3] == MPI_PROC_NULL)
    {
        new_grid[i_loc_max][1] = 0;
        new_grid[i_loc_max][j_loc_max] = 0;
    }
}

//! apply_neumann_bc
/*!
    Checks if neighbours are out of bounds,
    if so the physical boundary is set such that
    nabla u dot n  = 0.
*/
void algorithm::apply_neumann_bc()
{
    // LEFT
    if (neighbours[0] == MPI_PROC_NULL)
    {
        for (int i = 1; i <= i_loc_max; i++)
            new_grid[i][1] = new_grid[i][2];
    }
    // TOP
    if (neighbours[1] == MPI_PROC_NULL)
    {
        for (int j = 1; j <= j_loc_max; j++)
            new_grid[1][j] = new_grid[2][j];
    }
    // RIGHT
    if (neighbours[2] == MPI_PROC_NULL)
    {
        for (int i = 1; i <= i_loc_max; i++)
            new_grid[i][j_loc_w_ghost_max - 2] = new_grid[i][j_loc_w_ghost_max - 3];
    }
    // BOTTOM
    if (neighbours[3] == MPI_PROC_NULL)
    {
        for (int j = 1; j <= j_loc_max; j++)
            new_grid[i_loc_w_ghost_max - 2][j] = new_grid[i_loc_w_ghost_max - 3][j];
    }
}

//!
/*!
    Sets top row ghost cells with bottom (not ghost cells) rows 
    and vice versa. Only used for periodic boundaries and when 
    a prime number of processors are used.
*/
void algorithm::set_periodic_without_comms()
{
    for (int j = 1; j <= j_loc_max; j++)
    {
        grid[0][j] = grid[i_loc_w_ghost_max - 2][j];
        grid[i_loc_w_ghost_max - 1][j] = grid[1][j];
    }
}
//! grid_to_file
/*!
    grid_to_file takes a parameter out
    which is attached to the file name.
    It prints the grid to a file.
*/
void algorithm::grid_to_file(int out)
{
    stringstream fname;
    fstream f1;
    fname << "./out/output"
          << "_" << out << "_"
          << "id" << id << ".dat";
    f1.open(fname.str().c_str(), ios_base::out);
    f1 << m << "\t" << n << "\n";
    f1 << i_dom * j_loc_max << "\t" << j_dom * i_loc_max << "\n";
    for (int i = 1; i <= i_loc_max; i++)
    {
        for (int j = 1; j <= j_loc_max; j++)
            f1 << grid[i][j] << "\t";
        f1 << "\n";
    }
    f1.close();
}

//! run_dirchlet
/*!
    Calls dirchlet boundary specific functions and general simulaiton functions.
    The simulation is ran until t < tmax criterion is met.
*/
void algorithm::run_dirchlet()
{
    calc_neighbours("Dirchlet");
    setup_continuous_array(grid, grid1, i_loc_w_ghost_max, j_loc_w_ghost_max);
    setup_continuous_array(new_grid, new_grid1, i_loc_w_ghost_max, j_loc_w_ghost_max);
    setup_continuous_array(old_grid, old_grid1, i_loc_w_ghost_max, j_loc_w_ghost_max);
    create_datatypes(grid, grid1);
    direction();
    initial_condition();
    apply_dirchlet_bc();
    grid_to_file(out_cnt);
    out_cnt++;
    t_out += dt_out;
    while (t < t_max)
    {

        do_comms();
        do_iteration_ext_periodic_non_periodic_no_edge();
        apply_dirchlet_bc_loop();
        // swapping after all iterations on a grid are complete.
        std::swap(this->old_grid1, this->new_grid1);
        std::swap(this->old_grid, this->new_grid);
        std::swap(this->old_grid1, this->grid1);        //2D SWAP
        std::swap(this->old_grid, this->grid); 

        if (t_out <= t)
        {
            grid_to_file(out_cnt);
            out_cnt++;
            t_out += dt_out;
        }
        t += dt;
        it++;
    }
}
//! run_periodic
/*!
    Calls periodic boundary specific functions  and general simulaiton functions.
    The simulation is ran until t < tmax criterion is met.
*/
void algorithm::run_periodic()
{
    calc_neighbours("Periodic");
    setup_continuous_array(grid, grid1, i_loc_w_ghost_max, j_loc_w_ghost_max);
    setup_continuous_array(new_grid, new_grid1, i_loc_w_ghost_max, j_loc_w_ghost_max);
    setup_continuous_array(old_grid, old_grid1, i_loc_w_ghost_max, j_loc_w_ghost_max);
    create_datatypes(grid, grid1);
    direction();
    initial_condition();
    grid_to_file(out_cnt);
    out_cnt++;
    t_out += dt_out;
    while (t < t_max)
    {
        do_comms();
        ext_iteration_periodic();
        std::swap(this->old_grid1, this->new_grid1);
        std::swap(this->old_grid, this->new_grid);
        std::swap(this->old_grid1, this->grid1);        //2D SWAP
        std::swap(this->old_grid, this->grid); 
        if (t_out <= t)
        {
          grid_to_file(out_cnt);
          out_cnt++;
          t_out += dt_out;
        }
        t += dt;
        it++;
    }
}
//!run_neumann
/*!
    Calls neumann boundary specific functions  and general simulaiton functions.
    The simulation is ran until t < tmax criterion is met.
*/
void algorithm::run_neumann()
{
    calc_neighbours("Neumann");
    setup_continuous_array(grid, grid1, i_loc_w_ghost_max, j_loc_w_ghost_max);
    setup_continuous_array(new_grid, new_grid1, i_loc_w_ghost_max, j_loc_w_ghost_max);
    setup_continuous_array(old_grid, old_grid1, i_loc_w_ghost_max, j_loc_w_ghost_max);
    create_datatypes(grid, grid1);
    direction();
    initial_condition();
    grid_to_file(out_cnt);
    out_cnt++;
    t_out += dt_out;
    while (t < t_max)
    {
        do_comms();
        do_iteration_ext_periodic_non_periodic_no_edge();
        apply_neumann_bc();
	std::swap(this->old_grid1, this->new_grid1);
        std::swap(this->old_grid, this->new_grid);
        std::swap(this->old_grid1, this->grid1);	//2D SWAP
        std::swap(this->old_grid, this->grid);     //2D SWAP
        if (t_out <= t)
        {
            grid_to_file(out_cnt);
            out_cnt++;
            t_out += dt_out;
        }
        t += dt;
        it++;
    }
}

//!
/*!
    Calls periodic boundary and m = 1 specific functions  and general 
    simulaiton functions. The simulation is ran until t < tmax criterion is met.
*/
void algorithm::run_periodic_prime_processor()
{
    calc_neighbours("Periodic");
    setup_continuous_array(grid, grid1, i_loc_w_ghost_max, j_loc_w_ghost_max);
    setup_continuous_array(new_grid, new_grid1, i_loc_w_ghost_max, j_loc_w_ghost_max);
    setup_continuous_array(old_grid, old_grid1, i_loc_w_ghost_max, j_loc_w_ghost_max);
    create_datatypes(grid, grid1);
    direction();
    grid_to_file(out_cnt);
    out_cnt++;
   t_out += dt_out;
    initial_condition();
    while (t < t_max)
    {
        do_comms();
        set_periodic_without_comms();
        ext_iteration_periodic();
	std::swap(this->old_grid1, this->new_grid1);
        std::swap(this->old_grid, this->new_grid);
	std::swap(this->old_grid1, this->grid1);	//2D SWAP
        std::swap(this->old_grid, this->grid);     //2D SWAP
       if (t_out <= t)
       {
          grid_to_file(out_cnt);
          out_cnt++;
           t_out += dt_out;
       }
        t += dt;
        it++;
    }
}
