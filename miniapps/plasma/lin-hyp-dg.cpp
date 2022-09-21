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

/*
need to write class for DG Solver
*/

//////// START MAIN /////////
int main(int argc, char *argv[])
{

    // Initialise MPI and HYPRE
    Mpi::Init(argc, argv);
    int num_procs = Mpi::WorldSize();
    int myid = Mpi::WorldRank();
    Hypre::Init();
    bool verbose = (myid == 0);

    // parse command line options
    const char *mesh_file = "~/Documents/git-repos/RF-SciDAC/mfem-analysis/1D_tests/linear-hyp/slab_128.mesh";
    int precision 16;
    cout.precision(precision);

    OptionsParser args(argc,argv);
    args.AddOption(&mesh_file, "-m", "--mesh",
                    "Mesh file to use");
    args.Parse();
    if (!args.Good())
    {
        if (Mpi::Root)
        {
            args.PrintUsage(cout);
        }
        return 1;
    }
    if (Mpi::Root())
    {
        args.PrintOptions(cout);
    } 

/*
define mesh file: can be simple uniform grid to start e.g. slab_128.mesh.
*/

    // Read mesh file
    // Not including refinement for now but 
    Mesh mesh(mesh_file, 1, 1);
    const int dim = mesh.Dimension(); 
    /*
    What's the difference between defining the mesh dimension as 
        const int dim =  mesh.Dimension()
    and
        int dim = mesh->Dimension();
    ??
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
