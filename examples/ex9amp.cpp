//                       MFEM Example 9 - Parallel Version
//
// Compile with: make ex9p
//
// Sample runs:
//    mpirun -np 4 ex9p -m ../data/periodic-segment.mesh -p 0 -dt 0.005
//    mpirun -np 4 ex9p -m ../data/periodic-square.mesh -p 0 -dt 0.01
//    mpirun -np 4 ex9p -m ../data/periodic-hexagon.mesh -p 0 -dt 0.01
//    mpirun -np 4 ex9p -m ../data/periodic-square.mesh -p 1 -dt 0.005 -tf 9
//    mpirun -np 4 ex9p -m ../data/periodic-hexagon.mesh -p 1 -dt 0.005 -tf 9
//    mpirun -np 4 ex9p -m ../data/amr-quad.mesh -p 1 -rp 1 -dt 0.002 -tf 9
//    mpirun -np 4 ex9p -m ../data/amr-quad.mesh -p 1 -rp 1 -dt 0.02 -s 13 -tf 9
//    mpirun -np 4 ex9p -m ../data/star-q3.mesh -p 1 -rp 1 -dt 0.004 -tf 9
//    mpirun -np 4 ex9p -m ../data/star-mixed.mesh -p 1 -rp 1 -dt 0.004 -tf 9
//    mpirun -np 4 ex9p -m ../data/disc-nurbs.mesh -p 1 -rp 1 -dt 0.005 -tf 9
//    mpirun -np 4 ex9p -m ../data/disc-nurbs.mesh -p 2 -rp 1 -dt 0.005 -tf 9
//    mpirun -np 4 ex9p -m ../data/periodic-square.mesh -p 3 -rp 2 -dt 0.0025 -tf 9 -vs 20
//    mpirun -np 4 ex9p -m ../data/periodic-cube.mesh -p 0 -o 2 -rp 1 -dt 0.01 -tf 8
//    mpirun -np 3 ex9p -m ../data/amr-hex.mesh -p 1 -rs 1 -rp 0 -dt 0.005 -tf 0.5
//
// Device sample runs:
//    mpirun -np 4 ex9p -pa
//    mpirun -np 4 ex9p -ea
//    mpirun -np 4 ex9p -fa
//    mpirun -np 4 ex9p -pa -m ../data/periodic-cube.mesh
//    mpirun -np 4 ex9p -pa -m ../data/periodic-cube.mesh -d cuda
//    mpirun -np 4 ex9p -ea -m ../data/periodic-cube.mesh -d cuda
//    mpirun -np 4 ex9p -fa -m ../data/periodic-cube.mesh -d cuda
//
// Description:  This example code solves the time-dependent advection equation
//               du/dt + v.grad(u) = 0, where v is a given fluid velocity, and
//               u0(x)=u(0,x) is a given initial condition.
//
//               The example demonstrates the use of Discontinuous Galerkin (DG)
//               bilinear forms in MFEM (face integrators), the use of implicit
//               and explicit ODE time integrators, the definition of periodic
//               boundary conditions through periodic meshes, as well as the use
//               of GLVis for persistent visualization of a time-evolving
//               solution. Saving of time-dependent data files for visualization
//               with VisIt (visit.llnl.gov) and ParaView (paraview.org), as
//               well as the optional saving with ADIOS2 (adios2.readthedocs.io)
//               are also illustrated.
//
// mpirun -np 1 ./ex9amp -m ../miniapps/meshing/annulus-quad-o3.mesh -p 5 -s 23 -tf 8 -vs 10 -rs 0 -rp 0 -no-vis -o 3 -dt 0.025 -s0e .1 -t 1 -pt 1 -s1e -1e-3
//
#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

// Choice for the problem setup. The fluid velocity, initial condition and
// inflow boundary condition are chosen based on this parameter.
int problem;

// Velocity coefficient
void velocity_function(const Vector &x, Vector &v);

// Initial condition
double u0_function(const Vector &x, double t);

// Inflow boundary condition
double inflow_function(const Vector &x);

// Mesh bounding box
Vector bb_min, bb_max;

// Type of preconditioner for implicit time integrator
enum class PrecType : int
{
   ILU = 0,
   AIR = 1
};

#if MFEM_HYPRE_VERSION >= 21800
// Algebraic multigrid preconditioner for advective problems based on
// approximate ideal restriction (AIR). Most effective when matrix is
// first scaled by DG block inverse, and AIR applied to scaled matrix.
// See https://doi.org/10.1137/17M1144350.
class AIR_prec : public Solver
{
private:
   const HypreParMatrix *A;
   // Copy of A scaled by block-diagonal inverse
   HypreParMatrix A_s;

   HypreBoomerAMG *AIR_solver;
   int blocksize;

public:
   AIR_prec(int blocksize_) : AIR_solver(NULL), blocksize(blocksize_) { }

   void SetOperator(const Operator &op)
   {
      width = op.Width();
      height = op.Height();

      A = dynamic_cast<const HypreParMatrix *>(&op);
      MFEM_VERIFY(A != NULL, "AIR_prec requires a HypreParMatrix.")

      // Scale A by block-diagonal inverse
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
      // Scale the rhs by block inverse and solve system
      HypreParVector z_s;
      BlockInverseScale(A, NULL, &x, &z_s, blocksize,
                        BlockInverseScaleJob::RHS_ONLY);
      AIR_solver->Mult(z_s, y);
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
   mutable int min_its, max_its, cum_its, num_mult;
  
public:
   DG_Solver(HypreParMatrix &M_, HypreParMatrix &K_,
	     const FiniteElementSpace &fes,
             PrecType prec_type)
      : M(M_),
        K(K_),
        A(NULL),
        linear_solver(M.GetComm()),
        dt(-1.0),
	min_its(100),
	max_its(0),
	cum_its(0),
	num_mult(0)
   {
      int block_size = fes.GetFE(0)->GetDof();
      if (prec_type == PrecType::ILU)
      {
         prec = new BlockILU(block_size,
                             BlockILU::Reordering::MINIMUM_DISCARDED_FILL);
      }
      else if (prec_type == PrecType::AIR)
      {
#if MFEM_HYPRE_VERSION >= 21800
         prec = new AIR_prec(block_size);
#else
         MFEM_ABORT("Must have MFEM_HYPRE_VERSION >= 21800 to use AIR.\n");
#endif
      }
      linear_solver.iterative_mode = false;
      linear_solver.SetRelTol(1e-9);
      linear_solver.SetAbsTol(0.0);
      linear_solver.SetMaxIter(100);
      linear_solver.SetPrintLevel(0);
      linear_solver.SetPreconditioner(*prec);

      M.GetDiag(M_diag);
   }

   double GetMaximumTimeStep() const;
  
   void SetTimeStep(double dt_)
   {
      if (dt_ != dt)
      {
         dt = dt_;
         // Form operator A = M - dt*K
         delete A;
         A = Add(-dt, K, 0.0, K);
         SparseMatrix A_diag;
         A->GetDiag(A_diag);
         A_diag.Add(1.0, M_diag);
         // this will also call SetOperator on the preconditioner
         linear_solver.SetOperator(*A);
      }
   }

   void SetOperator(const Operator &op)
   {
      linear_solver.SetOperator(op);
   }

   virtual void Mult(const Vector &x, Vector &y) const
   {
      linear_solver.Mult(x, y);

      int num_its = linear_solver.GetNumIterations();

      min_its = std::min(min_its, num_its);
      max_its = std::max(max_its, num_its);
      cum_its += num_its;
      num_mult++;
   }
  
   ~DG_Solver()
   {
      delete prec;
      delete A;
   }

  void PrintIterationStats() const
  {
    mfem::out << "Number of iterations (min, max, avg): "
	      << min_its << ", " << max_its << ", " << cum_its / num_mult
	      << endl;
  }
};


/** A time-dependent operator for the right-hand side of the ODE. The DG weak
    form of du/dt = -v.grad(u) is M du/dt = K u + b, where M and K are the mass
    and advection matrices, and b describes the flow on the boundary. This can
    be written as a general ODE, du/dt = M^{-1} (K u + b), and this class is
    used to evaluate the right-hand side. */
class FE_Evolution : public TimeDependentOperator
{
private:
   OperatorHandle M, K;
   const Vector &b;
   Solver *M_prec;
   CGSolver M_solver;
   DG_Solver *dg_solver;

   mutable Vector z;

public:
   FE_Evolution(ParBilinearForm &_M, ParBilinearForm &_K, const Vector &_b,
                PrecType prec_type);

   virtual void Mult(const Vector &x, Vector &y) const;
   virtual void ImplicitSolve(const double dt, const Vector &x, Vector &k);

   virtual ~FE_Evolution();

  void PrintSolverStats() const
  {
    if (dg_solver)
      {
	dg_solver->PrintIterationStats();
      }
  }
};

class DGAdvDiffIntegrator : public BilinearFormIntegrator
{
protected:
  //double alpha;
  double s0_e, s1_e;
  Coefficient *tau;
  VectorCoefficient *beta;

  Vector shape1, shape2, ndshape1, ndshape2;
  DenseMatrix dshape1, dshape2;
  
public:
  DGAdvDiffIntegrator(VectorCoefficient & b, Coefficient &t,
		      double s0, double s1)
    :
    // alpha(a),
    s0_e(s0),
    s1_e(s1),
    tau(&t),
    beta(&b)
  { }
   using BilinearFormIntegrator::AssembleFaceMatrix;
   virtual void AssembleFaceMatrix(const FiniteElement &el1,
                                   const FiniteElement &el2,
                                   FaceElementTransformations &Trans,
                                   DenseMatrix &elmat)
  {
   int dim, ndof1, ndof2;

   double w, t1, t2, b1n, b2n, alpha1, alpha2;

   dim = el1.GetDim();
   Vector nor(dim);
   Vector B1(dim), B2(dim);
   
   MFEM_VERIFY(Trans.Elem2No >= 0, "Designed for interior faces")

   ndof1 = el1.GetDof();
   ndof2 = el2.GetDof();

   shape1.SetSize(ndof1);
   shape2.SetSize(ndof2);
   dshape1.SetSize(ndof1, dim);
   dshape2.SetSize(ndof2, dim);
   ndshape1.SetSize(ndof1);
   ndshape2.SetSize(ndof2);
   elmat.SetSize(ndof1 + ndof2);
   elmat = 0.0;

   Vector c1(dim), c2(dim), x(dim);
   Vector v12(dim), c1x(dim), xc2(dim);
   {
     IntegrationPoint cent;
     cent.x = 0.5; cent.y = 0.5;
     Trans.Elem1->SetIntPoint(&cent);
     Trans.Elem2->SetIntPoint(&cent);
     Trans.Elem1->Transform(cent, c1);
     Trans.Elem2->Transform(cent, c2);
     subtract(c2, c1, v12);
   }
   
   const IntegrationRule *ir = IntRule;
   if (ir == NULL)
   {
      int order;
      // Assuming order(u)==order(mesh)
      if (Trans.Elem2No >= 0)
         order = (min(Trans.Elem1->OrderW(), Trans.Elem2->OrderW()) +
                  2*max(el1.GetOrder(), el2.GetOrder()));
      else
      {
         order = Trans.Elem1->OrderW() + 2*el1.GetOrder();
      }
      if (el1.Space() == FunctionSpace::Pk)
      {
         order++;
      }
      ir = &IntRules.Get(Trans.GetGeometryType(), order);
      if (Trans.Elem2No < 0 && false)
	{
	  mfem::out << "DGTrace order " << order
		    << ", num pts = " << ir->GetNPoints() << std::endl; 
	}
   }

   for (int p = 0; p < ir->GetNPoints(); p++)
   {
      const IntegrationPoint &ip = ir->IntPoint(p);

      // Set the integration point in the face and the neighboring elements
      Trans.SetAllIntPoints(&ip);

      {
	Trans.Transform(ip, x);
      }
      
      // Access the neighboring elements' integration points
      // Note: eip2 will only contain valid data if Elem2 exists
      const IntegrationPoint &eip1 = Trans.GetElement1IntPoint();
      const IntegrationPoint &eip2 = Trans.GetElement2IntPoint();

      el1.CalcShape(eip1, shape1);
      el2.CalcShape(eip2, shape2);

      // el1.CalcDShape(eip1, dshape1);
      // el2.CalcDShape(eip2, dshape2);
      el1.CalcPhysDShape(*Trans.Elem1, dshape1);
      el2.CalcPhysDShape(*Trans.Elem1, dshape2);

      beta->Eval(B1, *Trans.Elem1, eip1);
      beta->Eval(B2, *Trans.Elem2, eip2);

      t1 = tau->Eval(*Trans.Elem1, eip1);
      t2 = tau->Eval(*Trans.Elem2, eip2);
      
      if (dim == 1)
      {
         nor(0) = 2*eip1.x - 1.0;
      }
      else
      {
         CalcOrtho(Trans.Jacobian(), nor);
	 // mfem::out << "|nor| " << nor.Norml2() << ", |J| " << Trans.Jacobian().Weight() << endl;
	 subtract(c2, x, xc2);
	 subtract(x, c1, c1x);
	 if ( v12 * nor <= 0.0)
	   {
	 mfem::out << "c1 (" << c1[0] << "," << c1[1] << ")"
		   << ", c2 (" << c2[0] << "," << c2[1] << ")"
		   << ", x (" << x[0] << "," << x[1] << ")"
		   << ", n (" << nor[0] << "," << nor[1] << ")"
		   << ", " << v12 * nor << ", " << c1x * nor << ", " << xc2 * nor << endl;
	   }
      }

      double detJ = nor.Norml2();
      
      b1n = B1 * nor;
      b2n = B2 * nor;

      alpha1 = 0.5 * (1.0 + copysign(t1, b1n));
      alpha2 = 1.0 - alpha1;

      dshape1.Mult(nor, ndshape1);
      dshape2.Mult(nor, ndshape2);
      
      // un = vu * nor;
      // a = 0.5 * alpha * un;
      // b = beta * fabs(un);
      // note: if |alpha/2|==|beta| then |a|==|b|, i.e. (a==b) or (a==-b)
      //       and therefore two blocks in the element matrix contribution
      //       (from the current quadrature point) are 0
      /*
      if (Trans.Elem2No < 0 && false)
	{
	  mfem::out << "DGTrace v = (" << vu[0] << ", " << vu[1] << "), n "
		    << nor[0] << ", " << nor[1] << ")" << std::endl;
	  mfem::out << "DGTrace " << un << ", a = " << a << ", b = " << b << std::endl;
	}
      */
      /*
      if (rho)
      {
         double rho_p;
         if (un >= 0.0 && ndof2)
         {
            rho_p = rho->Eval(*Trans.Elem2, eip2);
         }
         else
         {
            rho_p = rho->Eval(*Trans.Elem1, eip1);
         }
         a *= rho_p;
         b *= rho_p;
      }
      */
      w = (detJ * s0_e + (2.0 * alpha1 - 0.5) * b1n) * ip.weight;
      if (w != 0.0)
      {
         for (int i = 0; i < ndof1; i++)
            for (int j = 0; j < ndof1; j++)
            {
               elmat(i, j) += w * shape1(i) * shape1(j);
            }
      }

      w = (detJ * s0_e + alpha1 * b1n - (alpha2 - 0.5) * b2n) * ip.weight;
      if (w != 0.0)
      {
	for (int i = 0; i < ndof2; i++)
	  for (int j = 0; j < ndof1; j++)
	  {
	    elmat(ndof1+i, j) -= w * shape2(i) * shape1(j);
	  }
      }

      w = (detJ * s0_e + (alpha1 - 0.5) * b1n - alpha2 * b2n) * ip.weight;
      if (w != 0.0)
      {
	for (int i = 0; i < ndof1; i++)
	  for (int j = 0; j < ndof2; j++)
	  {
	    elmat(i, ndof1+j) -= w * shape1(i) * shape2(j);
	  }
      }
      
      w = (detJ * s0_e - (2.0 * alpha2 - 0.5) * b2n) * ip.weight;
      if (w != 0.0)
      {
	for (int i = 0; i < ndof2; i++)
	  for (int j = 0; j < ndof2; j++)
	  {
	    elmat(ndof1+i, ndof1+j) += w * shape2(i) * shape2(j);
	  }
      }
      /*
      if (s1_e != 0.0)
      {
	alpha1 = 0.5;
	alpha2 = 0.5;
	w = s1_e * ip.weight;
	for (int i = 0; i < ndof1; i++)
	  for (int j = 0; j < ndof1; j++)
          {
	    elmat(i, j) -= alpha1 * w * shape1(i) * ndshape1(j);
	    elmat(i, j) -= alpha1 * w * ndshape1(i) * shape1(j);
	  }

	for (int i = 0; i < ndof2; i++)
	  for (int j = 0; j < ndof1; j++)
	  {
	    elmat(ndof1+i, j) += alpha1 * w * shape2(i) * ndshape1(j);
	    elmat(ndof1+i, j) -= alpha2 * w * ndshape2(i) * shape1(j);
	  }

	for (int i = 0; i < ndof1; i++)
	  for (int j = 0; j < ndof2; j++)
	  {
	    elmat(i, ndof1+j) -= alpha2 * w * shape1(i) * ndshape2(j);
	    elmat(i, ndof1+j) += alpha1 * w * ndshape1(i) * shape2(j);
	  }
      
	for (int i = 0; i < ndof2; i++)
	  for (int j = 0; j < ndof2; j++)
	  {
	    elmat(ndof1+i, ndof1+j) += alpha2 * w * shape2(i) * ndshape2(j);
	    elmat(ndof1+i, ndof1+j) += alpha2 * w * ndshape2(i) * shape2(j);
	  }
      }
      */
      if (s1_e != 0.0)
      {
	w = s1_e * ip.weight / detJ;
	for (int i = 0; i < ndof1; i++)
	  for (int j = 0; j < ndof1; j++)
          {
	    elmat(i, j) += w * ndshape1(i) * ndshape1(j);
	  }

	for (int i = 0; i < ndof2; i++)
	  for (int j = 0; j < ndof1; j++)
	  {
	    elmat(ndof1+i, j) -= w * ndshape2(i) * ndshape1(j);
	  }

	for (int i = 0; i < ndof1; i++)
	  for (int j = 0; j < ndof2; j++)
	  {
	    elmat(i, ndof1+j) -= w * ndshape1(i) * ndshape2(j);
	  }
      
	for (int i = 0; i < ndof2; i++)
	  for (int j = 0; j < ndof2; j++)
	  {
	    elmat(ndof1+i, ndof1+j) += w * ndshape2(i) * ndshape2(j);
	  }
      }
   }
   elmat *= -1.0;
  }
};

int main(int argc, char *argv[])
{
   // 1. Initialize MPI.
   MPI_Session mpi;
   int num_procs = mpi.WorldSize();
   int myid = mpi.WorldRank();

   // 2. Parse command-line options.
   problem = 0;
   const char *mesh_file = "../data/periodic-hexagon.mesh";
   int ser_ref_levels = 2;
   int par_ref_levels = 0;
   int order = 3;
   bool pa = false;
   bool ea = false;
   bool fa = false;
   const char *device_config = "cpu";
   int ode_solver_type = 4;
   double t_final = 10.0;
   double dt = 0.01;
   double tau = 1.0;
   double s0_e = 1.0;
   double s1_e = 0.0;
   bool visualization = true;
   bool visit = false;
   bool paraview = false;
   bool adios2 = false;
   bool binary = false;
   int vis_steps = 5;
#if MFEM_HYPRE_VERSION >= 21800
   PrecType prec_type = PrecType::AIR;
#else
   PrecType prec_type = PrecType::ILU;
#endif
   int precision = 8;
   cout.precision(precision);

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&problem, "-p", "--problem",
                  "Problem setup to use. See options in velocity_function().");
   args.AddOption(&ser_ref_levels, "-rs", "--refine-serial",
                  "Number of times to refine the mesh uniformly in serial.");
   args.AddOption(&par_ref_levels, "-rp", "--refine-parallel",
                  "Number of times to refine the mesh uniformly in parallel.");
   args.AddOption(&order, "-o", "--order",
                  "Order (degree) of the finite elements.");
   args.AddOption(&pa, "-pa", "--partial-assembly", "-no-pa",
                  "--no-partial-assembly", "Enable Partial Assembly.");
   args.AddOption(&ea, "-ea", "--element-assembly", "-no-ea",
                  "--no-element-assembly", "Enable Element Assembly.");
   args.AddOption(&fa, "-fa", "--full-assembly", "-no-fa",
                  "--no-full-assembly", "Enable Full Assembly.");
   args.AddOption(&device_config, "-d", "--device",
                  "Device configuration string, see Device::Configure().");
   args.AddOption(&ode_solver_type, "-s", "--ode-solver",
                  "ODE solver: 1 - Forward Euler,\n\t"
                  "            2 - RK2 SSP, 3 - RK3 SSP, 4 - RK4, 6 - RK6,\n\t"
                  "            11 - Backward Euler,\n\t"
                  "            12 - SDIRK23 (L-stable), 13 - SDIRK33,\n\t"
                  "            22 - Implicit Midpoint Method,\n\t"
                  "            23 - SDIRK23 (A-stable), 24 - SDIRK34");
   args.AddOption(&t_final, "-tf", "--t-final",
                  "Final time; start time is 0.");
   args.AddOption(&dt, "-dt", "--time-step",
                  "Time step.");
   args.AddOption(&tau, "-t", "--tau",
                  "DG tau weight.");
   args.AddOption(&s0_e, "-s0e", "--s0-e",
                  "DG penalty.");
   args.AddOption(&s1_e, "-s1e", "--s1-e",
                  "DG penalty.");
   args.AddOption((int *)&prec_type, "-pt", "--prec-type", "Preconditioner for "
                  "implicit solves. 0 for ILU, 1 for pAIR-AMG.");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.AddOption(&visit, "-visit", "--visit-datafiles", "-no-visit",
                  "--no-visit-datafiles",
                  "Save data files for VisIt (visit.llnl.gov) visualization.");
   args.AddOption(&paraview, "-paraview", "--paraview-datafiles", "-no-paraview",
                  "--no-paraview-datafiles",
                  "Save data files for ParaView (paraview.org) visualization.");
   args.AddOption(&adios2, "-adios2", "--adios2-streams", "-no-adios2",
                  "--no-adios2-streams",
                  "Save data using adios2 streams.");
   args.AddOption(&binary, "-binary", "--binary-datafiles", "-ascii",
                  "--ascii-datafiles",
                  "Use binary (Sidre) or ascii format for VisIt data files.");
   args.AddOption(&vis_steps, "-vs", "--visualization-steps",
                  "Visualize every n-th timestep.");
   args.Parse();
   if (!args.Good())
   {
      if (mpi.Root())
      {
         args.PrintUsage(cout);
      }
      return 1;
   }
   if (mpi.Root())
   {
      args.PrintOptions(cout);
   }

   Device device(device_config);
   if (mpi.Root()) { device.Print(); }

   // 3. Read the serial mesh from the given mesh file on all processors. We can
   //    handle geometrically periodic meshes in this code.
   Mesh *mesh = new Mesh(mesh_file, 1, 1);
   int dim = mesh->Dimension();

   // 4. Define the ODE solver used for time integration. Several explicit
   //    Runge-Kutta methods are available.
   ODESolver *ode_solver = NULL;
   switch (ode_solver_type)
   {
      // Explicit methods
      case 1: ode_solver = new ForwardEulerSolver; break;
      case 2: ode_solver = new RK2Solver(1.0); break;
      case 3: ode_solver = new RK3SSPSolver; break;
      case 4: ode_solver = new RK4Solver; break;
      case 6: ode_solver = new RK6Solver; break;
      // Implicit (L-stable) methods
      case 11: ode_solver = new BackwardEulerSolver; break;
      case 12: ode_solver = new SDIRK23Solver(2); break;
      case 13: ode_solver = new SDIRK33Solver; break;
      // Implicit A-stable methods (not L-stable)
      case 22: ode_solver = new ImplicitMidpointSolver; break;
      case 23: ode_solver = new SDIRK23Solver; break;
      case 24: ode_solver = new SDIRK34Solver; break;
      case 534: ode_solver = new SDIRK534Solver; break;
      default:
         if (mpi.Root())
         {
            cout << "Unknown ODE solver type: " << ode_solver_type << '\n';
         }
         delete mesh;
         return 3;
   }

   // 5. Refine the mesh in serial to increase the resolution. In this example
   //    we do 'ser_ref_levels' of uniform refinement, where 'ser_ref_levels' is
   //    a command-line parameter. If the mesh is of NURBS type, we convert it
   //    to a (piecewise-polynomial) high-order mesh.
   for (int lev = 0; lev < ser_ref_levels; lev++)
   {
      mesh->UniformRefinement();
   }
   if (mesh->NURBSext)
   {
      mesh->SetCurvature(max(order, 1));
   }
   mesh->GetBoundingBox(bb_min, bb_max, max(order, 1));

   // 6. Define the parallel mesh by a partitioning of the serial mesh. Refine
   //    this mesh further in parallel to increase the resolution. Once the
   //    parallel mesh is defined, the serial mesh can be deleted.
   ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
   delete mesh;
   for (int lev = 0; lev < par_ref_levels; lev++)
   {
      pmesh->UniformRefinement();
   }

   // 7. Define the parallel discontinuous DG finite element space on the
   //    parallel refined mesh of the given polynomial order.
   DG_FECollection fec(order, dim, BasisType::GaussLobatto);
   ParFiniteElementSpace *fes = new ParFiniteElementSpace(pmesh, &fec);

   HYPRE_Int global_vSize = fes->GlobalTrueVSize();
   if (mpi.Root())
   {
      cout << "Number of unknowns: " << global_vSize << endl;
   }

   // 8. Set up and assemble the parallel bilinear and linear forms (and the
   //    parallel hypre matrices) corresponding to the DG discretization. The
   //    DGTraceIntegrator involves integrals over mesh interior faces.
   VectorFunctionCoefficient velocity(dim, velocity_function);
   FunctionCoefficient inflow(inflow_function);
   FunctionCoefficient u0(u0_function);
   ConstantCoefficient zero(0.0);
   ConstantCoefficient one(1.0);

   ParBilinearForm *m = new ParBilinearForm(fes);
   ParBilinearForm *k = new ParBilinearForm(fes);
   if (pa)
   {
      m->SetAssemblyLevel(AssemblyLevel::PARTIAL);
      k->SetAssemblyLevel(AssemblyLevel::PARTIAL);
   }
   else if (ea)
   {
      m->SetAssemblyLevel(AssemblyLevel::ELEMENT);
      k->SetAssemblyLevel(AssemblyLevel::ELEMENT);
   }
   else if (fa)
   {
      m->SetAssemblyLevel(AssemblyLevel::FULL);
      k->SetAssemblyLevel(AssemblyLevel::FULL);
   }

   m->AddDomainIntegrator(new MassIntegrator);
   // constexpr double alpha = -1.0;
   // constexpr double beta = -0.5;
   ConstantCoefficient tauCoef(tau);
   k->AddDomainIntegrator(new ConservativeConvectionIntegrator(velocity, -1.0));
   // k->AddDomainIntegrator(new MixedScalarWeakDivergenceIntegrator(velocity));
   k->AddInteriorFaceIntegrator(new DGAdvDiffIntegrator(velocity, tauCoef,
							s0_e, s1_e));
   k->AddBdrFaceIntegrator(new DGTraceIntegrator(velocity, 1.0, 0.5));

   ParLinearForm *b = new ParLinearForm(fes);
   // b->AddBdrFaceIntegrator(
   // new BoundaryFlowIntegrator(inflow, velocity, alpha));

   ParLinearForm *lf = new ParLinearForm(fes);
   lf ->AddDomainIntegrator(new DomainLFIntegrator(one));

   int skip_zeros = 0;
   m->Assemble();
   k->Assemble(skip_zeros);
   b->Assemble();
   lf->Assemble();
   m->Finalize();
   k->Finalize(skip_zeros);


   HypreParVector *B = b->ParallelAssemble();

   // 9. Define the initial conditions, save the corresponding grid function to
   //    a file and (optionally) save data in the VisIt format and initialize
   //    GLVis visualization.
   ParGridFunction *u = new ParGridFunction(fes);
   u->ProjectCoefficient(u0);
   HypreParVector *U = u->GetTrueDofs();
   double nrm0 = u->ComputeL2Error(zero);
   {
     double err = u->ComputeL2Error(u0);
     mfem::out << "Initial error " << err / nrm0 << endl;
   }

   double sum = (*lf)(*u);
   mfem::out << "Initial integral " << sum << endl;   

   {
      ostringstream mesh_name, sol_name;
      mesh_name << "ex9-mesh." << setfill('0') << setw(6) << myid;
      sol_name << "ex9-init." << setfill('0') << setw(6) << myid;
      ofstream omesh(mesh_name.str().c_str());
      omesh.precision(precision);
      pmesh->Print(omesh);
      ofstream osol(sol_name.str().c_str());
      osol.precision(precision);
      u->Save(osol);
   }

   // Create data collection for solution output: either VisItDataCollection for
   // ascii data files, or SidreDataCollection for binary data files.
   DataCollection *dc = NULL;
   if (visit)
   {
      if (binary)
      {
#ifdef MFEM_USE_SIDRE
         dc = new SidreDataCollection("Example9-Parallel", pmesh);
#else
         MFEM_ABORT("Must build with MFEM_USE_SIDRE=YES for binary output.");
#endif
      }
      else
      {
         dc = new VisItDataCollection("Example9-Parallel", pmesh);
         dc->SetPrecision(precision);
         // To save the mesh using MFEM's parallel mesh format:
         // dc->SetFormat(DataCollection::PARALLEL_FORMAT);
      }
      dc->RegisterField("solution", u);
      dc->SetCycle(0);
      dc->SetTime(0.0);
      dc->Save();
   }

   ParaViewDataCollection *pd = NULL;
   if (paraview)
   {
      pd = new ParaViewDataCollection("Example9P", pmesh);
      pd->SetPrefixPath("ParaView");
      pd->RegisterField("solution", u);
      pd->SetLevelsOfDetail(order);
      pd->SetDataFormat(VTKFormat::BINARY);
      pd->SetHighOrderOutput(true);
      pd->SetCycle(0);
      pd->SetTime(0.0);
      pd->Save();
   }

   // Optionally output a BP (binary pack) file using ADIOS2. This can be
   // visualized with the ParaView VTX reader.
#ifdef MFEM_USE_ADIOS2
   ADIOS2DataCollection *adios2_dc = NULL;
   if (adios2)
   {
      std::string postfix(mesh_file);
      postfix.erase(0, std::string("../data/").size() );
      postfix += "_o" + std::to_string(order);
      const std::string collection_name = "ex9-p-" + postfix + ".bp";

      adios2_dc = new ADIOS2DataCollection(MPI_COMM_WORLD, collection_name, pmesh);
      // output data substreams are half the number of mpi processes
      adios2_dc->SetParameter("SubStreams", std::to_string(num_procs/2) );
      // adios2_dc->SetLevelsOfDetail(2);
      adios2_dc->RegisterField("solution", u);
      adios2_dc->SetCycle(0);
      adios2_dc->SetTime(0.0);
      adios2_dc->Save();
   }
#endif

   socketstream sout;
   if (visualization)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      sout.open(vishost, visport);
      if (!sout)
      {
         if (mpi.Root())
            cout << "Unable to connect to GLVis server at "
                 << vishost << ':' << visport << endl;
         visualization = false;
         if (mpi.Root())
         {
            cout << "GLVis visualization disabled.\n";
         }
      }
      else
      {
         sout << "parallel " << num_procs << " " << myid << "\n";
         sout.precision(precision);
         sout << "solution\n" << *pmesh << *u;
         sout << "pause\n";
         sout << flush;
         if (mpi.Root())
            cout << "GLVis visualization paused."
                 << " Press space (in the GLVis window) to resume it.\n";
      }
   }

   ofstream ofsIts("ex9amp_its.dat");
   ofstream ofsErr("ex9amp_err.dat");
   
   // 10. Define the time-dependent evolution operator describing the ODE
   //     right-hand side, and perform time-integration (looping over the time
   //     iterations, ti, with a time-step dt).
   FE_Evolution adv(*m, *k, *B, prec_type);

   double t = 0.0;
   adv.SetTime(t);
   ode_solver->Init(adv);

   bool done = false;
   for (int ti = 0; !done; )
   {
      double dt_real = min(dt, t_final - t);
      ode_solver->Step(*U, t, dt_real);
      ti++;

      done = (t >= t_final - 1e-8*dt);

      if (done || ti % vis_steps == 0)
      {
         // 11. Extract the parallel grid function corresponding to the finite
         //     element approximation U (the local solution on each processor).
         *u = *U;

	 sum = (*lf)(*u);

	 u0.SetTime(t);
	 double err = u->ComputeL2Error(u0);

	 if (mpi.Root())
         {
            cout << "time step: " << ti << ", time: " << t
		 << ", integral " << sum
		 << ", range " << U->Min() << " to " << U->Max()
		 << ", error " << err / nrm0 << endl;
	    ofsErr << t << '\t' << err / nrm0 << endl;
         }

         if (visualization)
         {
            sout << "parallel " << num_procs << " " << myid << "\n";
            sout << "solution\n" << *pmesh << *u << flush;
         }

         if (visit)
         {
            dc->SetCycle(ti);
            dc->SetTime(t);
            dc->Save();
         }

         if (paraview)
         {
            pd->SetCycle(ti);
            pd->SetTime(t);
            pd->Save();
         }

#ifdef MFEM_USE_ADIOS2
         // transient solutions can be visualized with ParaView
         if (adios2)
         {
            adios2_dc->SetCycle(ti);
            adios2_dc->SetTime(t);
            adios2_dc->Save();
         }
#endif
      }
      /*
      if (t > ceil(t) - dt || t < floor(t) + dt)
      {
	u0.SetTime(t);
	double err = u->ComputeL2Error(u0);
	mfem::out << "Error at time " << t << " is " << err / nrm0 << endl;
      }
      */
   }

   ofsErr.close();

   adv.PrintSolverStats();
   
   // 12. Save the final solution in parallel. This output can be viewed later
   //     using GLVis: "glvis -np <np> -m ex9-mesh -g ex9-final".
   {
      *u = *U;
      ostringstream sol_name;
      sol_name << "ex9-final." << setfill('0') << setw(6) << myid;
      ofstream osol(sol_name.str().c_str());
      osol.precision(precision);
      u->Save(osol);
   }

   // 13. Free the used memory.
   delete U;
   delete u;
   delete B;
   delete b;
   delete k;
   delete m;
   delete fes;
   delete pmesh;
   delete ode_solver;
   delete pd;
#ifdef MFEM_USE_ADIOS2
   if (adios2)
   {
      delete adios2_dc;
   }
#endif
   delete dc;

   return 0;
}


// Implementation of class FE_Evolution
FE_Evolution::FE_Evolution(ParBilinearForm &_M, ParBilinearForm &_K,
                           const Vector &_b, PrecType prec_type)
   : TimeDependentOperator(_M.Height()), b(_b),
     M_solver(_M.ParFESpace()->GetComm()),
     z(_M.Height())
{
   if (_M.GetAssemblyLevel()==AssemblyLevel::LEGACYFULL)
   {
      M.Reset(_M.ParallelAssemble(), true);
      K.Reset(_K.ParallelAssemble(), true);
   }
   else
   {
      M.Reset(&_M, false);
      K.Reset(&_K, false);
   }

   M_solver.SetOperator(*M);

   Array<int> ess_tdof_list;
   if (_M.GetAssemblyLevel()==AssemblyLevel::LEGACYFULL)
   {
      HypreParMatrix &M_mat = *M.As<HypreParMatrix>();
      HypreParMatrix &K_mat = *K.As<HypreParMatrix>();
      HypreSmoother *hypre_prec = new HypreSmoother(M_mat, HypreSmoother::Jacobi);
      M_prec = hypre_prec;

      dg_solver = new DG_Solver(M_mat, K_mat, *_M.FESpace(), prec_type);
      double dtMax = dg_solver->GetMaximumTimeStep();
      mfem::out << "Max Stable Time Step: " << dtMax << std::endl;
   }
   else
   {
      M_prec = new OperatorJacobiSmoother(_M, ess_tdof_list);
      dg_solver = NULL;
   }

   M_solver.SetPreconditioner(*M_prec);
   M_solver.iterative_mode = false;
   M_solver.SetRelTol(1e-9);
   M_solver.SetAbsTol(0.0);
   M_solver.SetMaxIter(100);
   M_solver.SetPrintLevel(0);
}

// Solve the equation:
//    u_t = M^{-1}(Ku + b),
// by solving associated linear system
//    (M - dt*K) d = K*u + b
void FE_Evolution::ImplicitSolve(const double dt, const Vector &x, Vector &k)
{
   K->Mult(x, z);
   z += b;
   dg_solver->SetTimeStep(dt);
   dg_solver->Mult(z, k);
}

void FE_Evolution::Mult(const Vector &x, Vector &y) const
{
   // y = M^{-1} (K x + b)
   K->Mult(x, z);
   z += b;
   M_solver.Mult(z, y);
}

FE_Evolution::~FE_Evolution()
{
   delete M_prec;
   delete dg_solver;
}

double DG_Solver::GetMaximumTimeStep() const
{  
   HypreParVector * v0 = new HypreParVector(M);
   HypreParVector * v1 = new HypreParVector(M);
   HypreParVector * RHS = new HypreParVector(M);

   v0->Randomize(1234);

   int iter = 0, nstep = 20;
   double dt0 = 1.0, dt1 = 1.0, change = 1.0, ptol = 0.01;

   // Create Solver assuming no loss operators
   HypreDiagScale diag(M);
   HyprePCG pcg(M);
   pcg.SetTol(1e-12);
   pcg.SetMaxIter(200);
   pcg.SetPrintLevel(0);
   pcg.SetPreconditioner(diag);

   // Use power method to approximate the largest eigenvalue of the update
   // operator.
   while ( iter < nstep && change > ptol )
   {
      double normV0 = InnerProduct(*v0,*v0);
      *v0 /= sqrt(normV0);

      K.Mult(*v0,*RHS);
      pcg.Mult(*RHS,*v1);

      double lambda = InnerProduct(*v0,*v1);
      dt1 = 1.0/fabs(lambda);
      cout << "Norm V0 " << normV0 << ", lambda " << lambda << ", dt " << dt1 << endl;
      change = fabs((dt1-dt0)/dt0);
      dt0 = dt1;

      // if ( myid_ == 0 && logging_ > 1 )
      {
         cout << iter << ":  " << dt0 << " " << change << endl;
      }

      std::swap(v0, v1);

      iter++;
   }

   delete v0;
   delete v1;
   delete RHS;

   return dt0;
}

// Velocity coefficient
void velocity_function(const Vector &x, Vector &v)
{
   int dim = x.Size();

   // map to the reference [-1,1] domain
   Vector X(dim);
   for (int i = 0; i < dim; i++)
   {
      double center = (bb_min[i] + bb_max[i]) * 0.5;
      X(i) = 2 * (x(i) - center) / (bb_max[i] - bb_min[i]);
   }

   switch (problem)
   {
      case 0:
      {
         // Translations in 1D, 2D, and 3D
         switch (dim)
         {
            case 1: v(0) = 1.0; break;
            case 2: v(0) = sqrt(2./3.); v(1) = sqrt(1./3.); break;
            case 3: v(0) = sqrt(3./6.); v(1) = sqrt(2./6.); v(2) = sqrt(1./6.);
               break;
         }
         break;
      }
      case 1:
      case 2:
      {
         // Clockwise rotation in 2D around the origin
         const double w = M_PI/2;
         switch (dim)
         {
            case 1: v(0) = 1.0; break;
            case 2: v(0) = w*X(1); v(1) = -w*X(0); break;
            case 3: v(0) = w*X(1); v(1) = -w*X(0); v(2) = 0.0; break;
         }
         break;
      }
      case 3:
      {
         // Clockwise twisting rotation in 2D around the origin
         const double w = M_PI/2;
         double d = max((X(0)+1.)*(1.-X(0)),0.) * max((X(1)+1.)*(1.-X(1)),0.);
         d = d*d;
         switch (dim)
         {
            case 1: v(0) = 1.0; break;
            case 2: v(0) = d*w*X(1); v(1) = -d*w*X(0); break;
            case 3: v(0) = d*w*X(1); v(1) = -d*w*X(0); v(2) = 0.0; break;
         }
         break;
      }
      case 4:
      case 5:
      {
         // Clockwise rotation in 2D around the origin
         const double w = 2.0 * M_PI;
         switch (dim)
         {
            case 1: v(0) = 1.0; break;
            case 2: v(0) = -w*X(1); v(1) = w*X(0); break;
            case 3: v(0) = -w*X(1); v(1) = w*X(0); v(2) = 0.0; break;
         }
         break;
      }
   }
}

// Initial condition
double u0_function(const Vector &x, double t)
{
   int dim = x.Size();

   // map to the reference [-1,1] domain
   Vector X(dim);
   for (int i = 0; i < dim; i++)
   {
      double center = (bb_min[i] + bb_max[i]) * 0.5;
      X(i) = 2 * (x(i) - center) / (bb_max[i] - bb_min[i]);
   }

   switch (problem)
   {
      case 0:
      case 1:
      {
         switch (dim)
         {
            case 1:
               return exp(-40.*pow(X(0)-0.5,2));
            case 2:
            case 3:
            {
               double rx = 0.45, ry = 0.25, cx = 0., cy = -0.2, w = 10.;
               if (dim == 3)
               {
                  const double s = (1. + 0.25*cos(2*M_PI*X(2)));
                  rx *= s;
                  ry *= s;
               }
               return ( erfc(w*(X(0)-cx-rx))*erfc(-w*(X(0)-cx+rx)) *
                        erfc(w*(X(1)-cy-ry))*erfc(-w*(X(1)-cy+ry)) )/16;
            }
         }
      }
      case 2:
      {
         double x_ = X(0), y_ = X(1), rho, phi;
         rho = hypot(x_, y_);
         phi = atan2(y_, x_);
         return pow(sin(M_PI*rho),2)*sin(3*phi);
      }
      case 3:
      {
         const double f = M_PI;
         return sin(f*X(0))*sin(f*X(1));
      }
   case 4:
     {
       int n = 1;
       double a = 1e17;
       double b = 1e16;
       double phi = atan2(X[1], X[0]);
       return a + b * sin((double)n * phi - 2.0 * M_PI * t);
     }
   case 5:
     {
       int n = 1;
       double a = 1;
       double b = 0.1;
       double phi = atan2(X[1], X[0]);
       return a + b * sin((double)n * phi - 2.0 * M_PI * t);
     }
   }
   return 0.0;
}

// Inflow boundary condition (zero for the problems considered in this example)
double inflow_function(const Vector &x)
{
   switch (problem)
   {
      case 0:
      case 1:
      case 2:
      case 3: return 0.0;
   }
   return 0.0;
}