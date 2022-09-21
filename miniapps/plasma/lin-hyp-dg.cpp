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

enum class PrecType : int
{
    ILU = 0,
    AIR = 1
};

#if MFEM_HYPRE_VERSION >= 21800
class AIR_prec : public Solver
{
private:
    const HypreParMatrix *A;
    HypreParMatrix A_s;

    HypreBoomerAMG *AIR_solver;
    int blocksize;

public:
    AIR_prec(int blocksize_) : AIR_solver(NULL), blocksize(blocksize_) { }

    void SetOperator(const Operator &op)
    {
        width = op.width;
        height = op.height;

        A = dynamic_cast<const HyperParMatrix *>(&op);
        MFEM_VERIFY(A != NULL, "AIR_prec requires a HypreParMatrix.")

        BlockInverseScale(A, &A_s, NULL, NULL, blocksize,
                                BlockInverseScaleJob::MATRIX_ONLY);
        delete AIR_solver;
        AIR_solver = new HypreBoomerAMG(A_s);
        AIR_solver->SetAdvectiveOptions(1, "", "FA");
        AIR_solver->SetPrintLevel(0);
        AIR_solver->SetMaxLevels(50);
    }

    virtual void Mult(const Vector &x, Vector &y) const
    {
        HypreParVector z_s;
        BlockInverseScale(A, NULL, &x, &z_s, blocksize,
                            BlockInverseScaleJob::RHS_ONLY);
        AIR_solver->Mult(z_s,y);
    }

    ~AIR_prec()
    {
        delete AIR_solver;
    }
};
#endif

class DG_Solver : public Solver
{
private:
    HypreParMatrix &M, &K;
    SparseMatrix M_diag;
    HypreParMatrix *A;
    GMRESSolver linear_solver;
    Solver *prec;
    double dt;

public:
    DG_Solver(HypreParMatrix &M_, HypreParMatrix &K_, const FiniteElementSpace &fes,
                PrecType prec_type)
        : M(M_),
          K(K_),
          A(NULL),
          linear_solver(M.GetComm()),
          dt(-1.0)
    {

        int block_size = fes.GetFE(0)->GetDof();

#if MFEM_HYPRE_VERSION >= 21800
        prec = new AIR_prec(block_size);
#else
        MFEM_ABORT("Must have MFEM_HYPRE_VERSION >= 21800 to use AIR.\n");
#endif
        
        linear_solver.iterative_mode = false;
        linear_solver.SetRelTol(1e-9);
        linear_solver.SetAbsTol(0);
        linear_solver.SetMaxIter(100);
        linear_solver.SetPrintLevel(0);
        linear_solver.SetPreconditioner(*prec);
        
        M.GetDiag(M_diag);

    } 

    void SetTimeStep(double dt_)
    {
        if (dt_ !=dt)
        {
            dt = dt_;
            delete A;
            SparseMatrix A_diag;
            A->GetDiag(A_diag);
            A_diag.Add(1.0, M_diag);
            linear_solver.SetOperator(*A);

        }
    }

    void SetOperator(const Operator &op)
    {
        linear_solver.SetOperator(op);
    }

    virtual void Mult(const Vector &x, Vector &y) const
    {
        linear_solver.Mult(x,y);
    }

    ~DG_Solver()
    {
        delete prec;
        delete A;
    }
};

class FE_Evolution : public TimeDependentOperator
{
private:
    OperatorHandle M, K;
    const Vector &b;
    Solver *M_prec;
    CGSolver M_Solver;
    DG_Solver *dg_solver;

    mutable Vector z;

public:
    FE_Evolution(ParBilinearForm &M_, ParBilinearForm &K_, const Vector &b_);

    virtual void Mult(const Vector &x, Vector &y) const;
    virtual void ImplicitSolve(const double dt, const Vector &x, const Vector &k);

    virtual ~FE_Evolution();    
};

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
    int ser_ref_levels = 0;
    int par_ref_levels = 0;
#if MFEM_HYPRE_VERSION >= 21800
    PrecType prec_type = PrecType::AIR;
#endif
    int precision = 16;
    cout.precision(precision);

    OptionsParser args(argc,argv);
    args.AddOption(&mesh_file, "-m", "--mesh",
                    "Mesh file to use");
    args.AddOption(&ser_ref_levels, "-rs", "--refine-serial",
                    "Number of times to refine the mesh uniformly in serial.");
    args.AddOption(&par_ref_levels, "-rp", "--refine-parallel",
                    "Number of times to refine the mesh uniformly in parallel.");

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
    Mesh mesh(mesh_file, 1, 1);
    const int dim = mesh.Dimension(); 
    /*
    What's the difference between defining the mesh dimension as 
        const int dim =  mesh.Dimension()
    and
        int dim = mesh->Dimension();
    ??
    */
    
    // define ode solver: start w simple backward euler.
    ODESolver *ode_solver = new BackwardEulerSolver;

    // refine mesh in serial
    for (int lev = 0; lev < ser_ref_levels; lev++)
    {
        mesh.UniformRefinement();
    }

    // define parallel mesh.
    ParMesh pmesh(MPI_COMM_WORLD, mesh);
    mesh.Clear();
    for (int lev = 0; lev < par_ref_levels; lev++)
    {
        pmesh.UniformRefinement();
    }

        

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
