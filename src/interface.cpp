#include "interface.h"

//! Constructor
/*!
   Calls automatically gatherdetails function 
   upon declaration type of an object.
*/
interface::interface()
{
   gatherdetails();
}
//! Destructor
interface::~interface() {}
//! gatherdetails
/*!
   Opens a file, if it cannot be opened, an error message is 
   returned and the program exits. If opened successfully, 
   each even row number including line 0 is skipped due to 
   being a description of what variable values should be 
   inputted the line below. The values are read and assigned
   to the appropriate variable names.
*/
void interface::gatherdetails()
{
   std::string filename = "./input.txt";
   std::ifstream file(filename);

   if (!file)
   {
      std::cerr << "File could not be opened !\n";
      MPI_Finalize();
      exit(0);
   }

   std::string line;
   getline(file, line);
   file >> this->imax >> this->jmax;

   getline(file, line);
   file >> this->x_max >> this->y_max >> this->c;

   getline(file, line);
   file >> this->tmax >> this->dt_out >> this->t_out;

   getline(file, line);
   file >> this->boundary_condition;

   if (c < 0)
   {
      std::cerr << "c must be a real positive number";
      MPI_Finalize();
      exit(0);
   }
   double dx = this->x_max / ((double)this->imax - 1);
   double dy = this->y_max / ((double)this->jmax - 1);
   this->dt = 0.1 * std::min(dx, dy) / this->c;
   if ((dt + 2) < (dx / c))
   {
      std::cerr << "Numerical Scheme is unstable";
      MPI_Finalize();
      exit(0);
   }
}
