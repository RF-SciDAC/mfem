/*
Solve linear hyperbolic system: 

\partial_t n + \partial_x q = 0
\partial_t q + \partial_x n = 0

or diagonalised ODE:

\partial_t w^+ + \partial_x w^+ = 0
\partial_t w^- - \partial_x w^- = 0

*/

#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <algorithm>

using namespace std;
using namespace mfem;

/*
define any constants
*/

// advection speed 
const double vs = 1.0;

// define number of equations ? (taken from ex 14)
const int num_equations = 2;

//////// START MAIN /////////
int main()
{

    // Initialise MPI and HYPRE
    Mpi::Init
/*
define mesh file: can be simple uniform grid to start e.g. slab_128.mesh.
*/
    


/*
do we need mesh refinement? might be good to include for future...?
*/

/*
need to write class for DG Solver
*/

/*
write class for FE evolution

DG weak form for linear hyperbolic advection case is same as ex5, e.g. dn/dt = -dq/dx is M dn/dt = K q + b. 
For the coupled case, the block form would be...

        D = [M  K]
            [K  M]

and     u = [n]
            [q]
*/

}
