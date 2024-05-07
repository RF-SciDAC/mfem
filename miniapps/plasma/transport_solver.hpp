// Copyright (c) 2010-2023, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#ifndef MFEM_TRANSPORT_SOLVER
#define MFEM_TRANSPORT_SOLVER

#include "../common/fem_extras.hpp"
#include "../common/pfem_extras.hpp"
#include "../common/mesh_extras.hpp"
#include "plasma.hpp"
#include "flux_vectors.hpp"
#include "transport_coefs.hpp"

#ifdef MFEM_USE_MPI

namespace mfem
{

namespace plasma
{

namespace transport
{

/** Multispecies Electron-Ion Collision Time in seconds
 Te is the electron temperature in eV
 ns is the number of ion species
 ni is the density of ions (assuming ni=ne) in particles per meter^3
 zi is the charge number of the ion species
 lnLambda is the Coulomb Logarithm
*/
//double tau_e(double Te, int ns, double * ni, int * zi, double lnLambda);

/** Multispecies Ion-Ion Collision Time in seconds
   ma is the ion mass in a.m.u
   Ta is the ion temperature in eV
   ion is the ion species index for the desired collision time
   ns is the number of ion species
   ni is the density of ions (assuming ni=ne) in particles per meter^3
   zi is the charge number of the ion species
   lnLambda is the Coulomb Logarithm
*/
//double tau_i(double ma, double Ta, int ion, int ns, double * ni, int * zi,
//             double lnLambda);

/**
  Thermal diffusion coefficient along B field for electrons
  Return value is in m^2/s.
   Te is the electron temperature in eV
   ns is the number of ion species
   ni is the density of ions (assuming ni=ne) in particles per meter^3
   zi is the charge number of the ion species
*/
/*
inline double chi_e_para(double Te, int ns, double * ni, int * zi)
{
 // The factor of q_ is included to convert Te from eV to Joules
 return 3.16 * (q_ * Te / me_kg_) * tau_e(Te, ns, ni, zi, 17.0);
}
*/
/**
  Thermal diffusion coefficient perpendicular to B field for electrons
  Return value is in m^2/s.
*/
/*
inline double chi_e_perp()
{
 // The factor of q_ is included to convert Te from eV to Joules
 return 1.0;
}
*/
/**
  Thermal diffusion coefficient perpendicular to both B field and
  thermal gradient for electrons.
  Return value is in m^2/s.
   Te is the electron temperature in eV
   ni is the density of ions (assuming ni=ne) in particles per meter^3
   z is the charge number of the ion species
*/
/*
inline double chi_e_cross()
{
 // The factor of q_ is included to convert Te from eV to Joules
 return 0.0;
}
*/
/**
  Thermal diffusion coefficient along B field for ions
  Return value is in m^2/s.
   ma is the ion mass in a.m.u.
   Ta is the ion temperature in eV
   ion is the ion species index for the desired coefficient
   ns is the number of ion species
   nb is the density of ions in particles per meter^3
   zb is the charge number of the ion species
*/
/*
inline double chi_i_para(double ma, double Ta,
                       int ion, int ns, double * nb, int * zb)
{
 // The factor of q_ is included to convert Ta from eV to Joules
 // The factor of u_ is included to convert ma from a.m.u to kg
 return 3.9 * (q_ * Ta / (ma * kg_per_amu_ ) ) *
        tau_i(ma, Ta, ion, ns, nb, zb, 17.0);
}
*/
/**
  Thermal diffusion coefficient perpendicular to B field for ions
  Return value is in m^2/s.
*/
/*
inline double chi_i_perp()
{
 // The factor of q_ is included to convert Ti from eV to Joules
 // The factor of u_ is included to convert mi from a.m.u to kg
 return 1.0;
}
*/
/**
  Thermal diffusion coefficient perpendicular to both B field and
  thermal gradient for ions
  Return value is in m^2/s.
*/
/*
inline double chi_i_cross()
{
 // The factor of q_ is included to convert Ti from eV to Joules
 // The factor of u_ is included to convert mi from a.m.u to kg
 return 0.0;
}
*/
/**
  Viscosity coefficient along B field for electrons
  Return value is in (a.m.u)*m^2/s.
   ne is the density of electrons in particles per meter^3
   Te is the electron temperature in eV
   ns is the number of ion species
   ni is the density of ions (assuming ni=ne) in particles per meter^3
   zi is the charge number of the ion species
*/
/*
inline double eta_e_para(double ne, double Te, int ns, double * ni, int * zi)
{
 // The factor of q_ is included to convert Te from eV to Joules
 // The factor of u_ is included to convert from kg to a.m.u
 return 0.73 * ne * (q_ * Te / kg_per_amu_) * tau_e(Te, ns, ni, zi, 17.0);
}
*/
/**
  Viscosity coefficient along B field for ions
  Return value is in (a.m.u)*m^2/s.
   ma is the ion mass in a.m.u.
   Ta is the ion temperature in eV
   ion is the ion species index for the desired coefficient
   ns is the number of ion species
   nb is the density of ions in particles per meter^3
   zb is the charge number of the ion species
*/
/*
inline double eta_i_para(double ma, double Ta,
                       int ion, int ns, double * nb, int * zb)
{
 // The factor of q_ is included to convert Ti from eV to Joules
 // The factor of u_ is included to convert from kg to a.m.u
 return 0.96 * nb[ion] * (q_ * Ta / kg_per_amu_) *
        tau_i(ma, Ta, ion, ns, nb, zb, 17.0);
}
*/

void ElementOrder(ParFiniteElementSpace &fes, Vector &elemOrder);
void DiscontinuitySensor(GridFunction &u, double uRef, double alpha,
                         Vector &error);
void ParallelMeshSpacing(ParFiniteElementSpace &fes,
                         VectorCoefficient &B3, Vector &h);

struct CoefficientByAttr
{
   Array<int> attr;
   Coefficient * coef;
   bool ownCoef;
};

struct CoefficientsByAttr
{
   Array<int> attr;
   Array<Coefficient*> coefs;
   Array<bool> ownCoefs;
};

class AdvectionDiffusionBC
{
public:
   enum BCType {DIRICHLET_BC, NEUMANN_BC, ROBIN_BC, OUTFLOW_BC};

private:
   Array<CoefficientByAttr*>  dbc; // Dirichlet BC data
   Array<CoefficientByAttr*>  nbc; // Neumann BC data
   Array<CoefficientsByAttr*> rbc; // Robin BC data
   Array<CoefficientByAttr*>  obc; // Outflow BC data
   mutable Array<int>  hbc_attr; // Homogeneous Neumann BC boundary attributes
   Array<int>  dbc_attr; // Dirichlet BC boundary attributes
   Array<int>  obc_attr; // Outflow BC boundary attributes

   std::set<int> bc_attr;
   const Array<int> & bdr_attr;

   common::CoefFactory * coefFact;

   void ReadBCs(std::istream &input);

   void ReadAttr(std::istream &input,
                 BCType bctype,
                 Array<int> &attr);

   void ReadCoefByAttr(std::istream &input,
                       BCType bctype,
                       CoefficientByAttr &cba);

   void ReadCoefsByAttr(std::istream &input,
                        BCType bctype,
                        CoefficientsByAttr &cba);

public:
   AdvectionDiffusionBC(const Array<int> & bdr)
      : bdr_attr(bdr), coefFact(NULL) {}

   AdvectionDiffusionBC(const Array<int> & bdr,
                        common::CoefFactory &cf, std::istream &input)
      : bdr_attr(bdr), coefFact(&cf) { ReadBCs(input); }

   ~AdvectionDiffusionBC();

   void SetTime(double t) const;

   static const char * GetBCTypeName(BCType bctype);

   void LoadBCs(common::CoefFactory &cf, std::istream &input)
   { coefFact = &cf; ReadBCs(input); }

   // Enforce u = val on boundaries with attributes in bdr
   void AddDirichletBC(const Array<int> & bdr, Coefficient &val);

   // Enforce du/dn = val on boundaries with attributes in bdr
   void AddNeumannBC(const Array<int> & bdr, Coefficient &val);

   // Enforce du/dn + a u = b on boundaries with attributes in bdr
   void AddRobinBC(const Array<int> & bdr, Coefficient &a, Coefficient &b);

   // Allows restricted outflow of the fluid through the boundary
   /** An outflow boundary condition is zero on portions of the
       boundary where the advection is directed into the domain. On
       portions where the advection is directed outward a val = 1
       would allow all incident fluid to flow out of the domain. If
       val < 1 the outflow is restricted leading to a buildup of fluid
       at the boundary.
   */
   void AddOutflowBC(const Array<int> & bdr, Coefficient &val);

   const Array<CoefficientByAttr*> & GetDirichletBCs() const { return dbc; }
   const Array<CoefficientByAttr*> & GetNeumannBCs() const { return nbc; }
   const Array<CoefficientsByAttr*> & GetRobinBCs() const { return rbc; }
   const Array<CoefficientByAttr*> & GetOutflowBCs() const { return obc; }

   const Array<int> & GetHomogeneousNeumannBDR() const;
   const Array<int> & GetDirichletBDR() const { return dbc_attr; }
   const Array<int> & GetOutflowBDR() const { return obc_attr; }
};

/** A RecyclingBC describes recombination at a boundary

  In a Recycling boundary condition an ion species recombines with
  electrons contained within the surface of the domain boundary and
  the resulting neutral atoms are added to the population of neutrals.

  We will assume that diffusion into the wall can be neglected
  i.e. only advection of ions towards the wall will lead to
  recycling. It is possible that some fraction of ions will remain
  ionized. The `ion_frac` coefficient should return the fraction
  (between 0 and 1) of incident ions which will be absorbed by the
  boundary. The `neu_frac` coefficient should return the fraction of
  incident ions which will be recycled as neutrals. In summary
  `ion_frac` and `n eu_frac` should be chosen so that:
   0 <= neu_frac <= ion_frac <= 1.

*/
class RecyclingBC
{
private:
   int ion_index;
   int vel_index;
   int neu_index;

   Array<CoefficientsByAttr*> bc; // Recycling BC data

   common::CoefFactory * coefFact;

   // void ReadBCs(std::istream &input);
   void ReadBC(std::istream &input);

   void ReadAttr(std::istream &input,
                 Array<int> &attr);

   void ReadCoefsByAttr(std::istream &input,
                        CoefficientsByAttr &cba);

public:
   RecyclingBC()
      : coefFact(NULL) {}

   // RecyclingBC(common::CoefFactory &cf, std::istream &input)
   //  : coefFact(&cf) { ReadBCs(input); }

   ~RecyclingBC();

   void SetTime(double t) const;

   void LoadBCs(common::CoefFactory &cf, std::istream &input)
   { coefFact = &cf; ReadBC(input); }

   void AddRecyclingBC(int ion, int vel, int neu, const Array<int> & bdr,
                       Coefficient & ion_frac, Coefficient & neu_frac);

   int GetIonDensityIndex() const { return ion_index; }
   int GetIonVelocityIndex() const { return vel_index; }
   int GetNeutralDensityIndex() const { return neu_index; }

   const Array<CoefficientsByAttr*> & GetRecyclingBCs() const { return bc; }
};

class CoupledBCs
{
public:
   enum BCType {RECYCLING_BC};

private:
   Array<RecyclingBC*> rbcs_;

   void ReadBCs(common::CoefFactory &cf, std::istream &input);

public:
   CoupledBCs() {}
   ~CoupledBCs();

   void LoadBCs(common::CoefFactory &cf, std::istream &input)
   { ReadBCs(cf, input); }

   int GetNumRecyclingBCs() const { return rbcs_.Size(); }
   RecyclingBC & GetRecyclingBC(int i) { return *rbcs_[i]; }
   const RecyclingBC & GetRecyclingBC(int i) const { return *rbcs_[i]; }
};

class TransportBCs
{
private:
   int neqn_;
   Array<AdvectionDiffusionBC*> bcs_;
   const Array<int> bdr_attr_;

   CoupledBCs cbcs_;

   void ReadBCs(common::CoefFactory &cf, std::istream &input);

public:
   TransportBCs(const Array<int> & bdr_attr, int neqn);

   TransportBCs(const Array<int> & bdr_attr, int neqn,
                common::CoefFactory &cf, std::istream &input);

   ~TransportBCs();

   void LoadBCs(common::CoefFactory &cf, std::istream &input)
   { ReadBCs(cf, input); }

   AdvectionDiffusionBC & operator()(int i) { return *bcs_[i]; }
   const AdvectionDiffusionBC & operator()(int i) const { return *bcs_[i]; }

   AdvectionDiffusionBC & operator[](int i) { return *bcs_[i]; }
   const AdvectionDiffusionBC & operator[](int i) const { return *bcs_[i]; }

   CoupledBCs & GetCoupledBCs() { return cbcs_; }
   const CoupledBCs & GetCoupledBCs() const { return cbcs_; }
};
/*
class GeneralCoefficient
{
public:
  enum GenCoefType {SCALAR_COEF, VECTOR_COEF, MATRIX_COEF};

private:
   CoefFactory * coefFact;

   GenCoefType type;
   Coefficient       * sCoef;
   VectorCoefficient * vCoef;
   MatrixCoefficient * mCoef;

   void ReadCoef(std::istream &input);

public:
  GeneralCoefficient()
    : coefFact(NULL), sCoef(NULL), vCoef(NULL), mCoef(NULL) {}
  GeneralCoefficient(CoefFactory &cf, std::istream &input)
    : coefFact(&cf), sCoef(NULL), vCoef(NULL), mCoef(NULL)
  { ReadCoef(input); }

   void LoadCoef(CoefFactory &cf, std::istream &input)
   { coefFact = &cf; ReadCoef(input); }

   void AddCoefficient(Coefficient &c) { type = SCALAR_COEF; sCoef = &c; }
   void AddCoefficient(VectorCoefficient &c) { type = VECTOR_COEF; vCoef = &c; }
   void AddCoefficient(MatrixCoefficient &c) { type = MATRIX_COEF; mCoef = &c; }

   GenCoefType GetCoefficientType() const { return type; }
   Coefficient * GetCoefficient() const { return sCoef; }
   VectorCoefficient * GetVectorCoefficient() const { return vCoef; }
   MatrixCoefficient * GetMatrixCoefficient() const { return mCoef; }
};
*/
class TransportICs
{
private:
   int neqn_;
   Array<Coefficient *> ics_;
   Array<bool> own_ics_;

   void ReadICs(common::CoefFactory &cf, std::istream &input);

public:
   TransportICs(int neqn)
      : neqn_(neqn),
        ics_(neqn),
        own_ics_(neqn)
   {
      ics_ = NULL;
      own_ics_ = false;
   }

   TransportICs(int neqn, common::CoefFactory &cf, std::istream &input);

   ~TransportICs()
   {
      for (int i=0; i<neqn_; i++)
      {
         if (own_ics_[i]) { delete ics_[i]; }
      }
   }

   void LoadICs(common::CoefFactory &cf, std::istream &input)
   { ReadICs(cf, input); }

   void SetOwnership(int i, bool own) { own_ics_[i] = own; }

   Coefficient *& operator()(int i) { return ics_[i]; }
   const Coefficient * operator()(int i) const { return ics_[i]; }

   Coefficient *& operator[](int i) { return ics_[i]; }
   const Coefficient * operator[](int i) const { return ics_[i]; }
};

class TransportExactSolutions
{
private:
   int neqn_;
   Array<Coefficient *> ess_;
   Array<bool> own_ess_;

   void Read(common::CoefFactory &cf, std::istream &input);

public:
   TransportExactSolutions(int neqn)
      : neqn_(neqn),
        ess_(neqn),
        own_ess_(neqn)
   {
      ess_ = NULL;
      own_ess_ = false;
   }

   TransportExactSolutions(int neqn, common::CoefFactory &cf,
                           std::istream &input);

   ~TransportExactSolutions()
   {
      for (int i=0; i<neqn_; i++)
      {
         if (own_ess_[i]) { delete ess_[i]; }
      }
   }

   void LoadExactSolutions(common::CoefFactory &cf, std::istream &input)
   { Read(cf, input); }

   void SetOwnership(int i, bool own) { own_ess_[i] = own; }

   Coefficient *& operator()(int i) { return ess_[i]; }
   const Coefficient * operator()(int i) const { return ess_[i]; }

   Coefficient *& operator[](int i) { return ess_[i]; }
   const Coefficient * operator[](int i) const { return ess_[i]; }
};

class StateVariableCoef;
class StateVariableVecCoef;
class StateVariableMatCoef;

class EqnCoefficients
{
protected:

   Array<StateVariableCoef *> sCoefs_;
   Array<StateVariableVecCoef *> vCoefs_;
   Array<StateVariableMatCoef *> mCoefs_;

   std::vector<std::string> sCoefNames_;
   std::vector<std::string> vCoefNames_;
   std::vector<std::string> mCoefNames_;

   common::CoefFactory * coefFact;

   virtual void ReadCoefs(std::istream &input);

public:
   EqnCoefficients(int nSCoefs, int nVCoefs = 0, int nMCoefs = 0)
      : sCoefs_(nSCoefs), vCoefs_(nVCoefs), mCoefs_(nMCoefs),
        sCoefNames_(nSCoefs), vCoefNames_(nVCoefs), mCoefNames_(nMCoefs)
   {
      sCoefs_ = NULL;
      vCoefs_ = NULL;
      mCoefs_ = NULL;
   }

   virtual ~EqnCoefficients() {}

   void LoadCoefs(common::CoefFactory &cf, std::istream &input)
   { coefFact = &cf; ReadCoefs(input); }

   StateVariableCoef *& operator()(int i) { return sCoefs_[i]; }
   const StateVariableCoef * operator()(int i) const { return sCoefs_[i]; }

   StateVariableCoef *& operator[](int i) { return sCoefs_[i]; }
   const StateVariableCoef * operator[](int i) const { return sCoefs_[i]; }

   StateVariableCoef *& GetScalarCoefficient(int i) { return sCoefs_[i]; }
   const StateVariableCoef * GetScalarCoefficient(int i) const
   { return sCoefs_[i]; }

   StateVariableVecCoef *& GetVectorCoefficient(int i) { return vCoefs_[i]; }
   const StateVariableVecCoef * GetVectorCoefficient(int i) const
   { return vCoefs_[i]; }

   StateVariableMatCoef *& GetMatrixCoefficient(int i) { return mCoefs_[i]; }
   const StateVariableMatCoef * GetMatrixCoefficient(int i) const
   { return mCoefs_[i]; }
};

class NeutralDensityCoefs : public EqnCoefficients
{
public:
   enum sCoefNames {DIFFUSION_COEF = 0, SOURCE_COEF, NUM_SCALAR_COEFS};

   NeutralDensityCoefs();
};

class IonDensityCoefs : public EqnCoefficients
{
public:
   enum sCoefNames {PERP_DIFFUSION_COEF = 0, PARA_DIFFUSION_COEF,
                    SOURCE_COEF, NUM_SCALAR_COEFS
                   };

   IonDensityCoefs();
};

class IonMomentumCoefs : public EqnCoefficients
{
public:
   enum sCoefNames {PARA_DIFFUSION_COEF = 0, PERP_DIFFUSION_COEF,
                    SOURCE_COEF, NUM_SCALAR_COEFS
                   };

   IonMomentumCoefs();
};

class IonStaticPressureCoefs : public EqnCoefficients
{
public:
   enum sCoefNames {PARA_DIFFUSION_COEF = 0, PERP_DIFFUSION_COEF,
                    SOURCE_COEF, NUM_SCALAR_COEFS
                   };

   IonStaticPressureCoefs();
};

class ElectronStaticPressureCoefs : public EqnCoefficients
{
public:
   enum sCoefNames {PERP_DIFFUSION_COEF = 0, PARA_DIFFUSION_COEF,
                    SOURCE_COEF, NUM_SCALAR_COEFS
                   };

   ElectronStaticPressureCoefs();
};

class IonTotalEnergyCoefs : public EqnCoefficients
{
public:
   enum sCoefNames {PARA_DIFFUSION_COEF = 0, PERP_DIFFUSION_COEF,
                    SOURCE_COEF, NUM_SCALAR_COEFS
                   };

   IonTotalEnergyCoefs();
};

class ElectronTotalEnergyCoefs : public EqnCoefficients
{
public:
   enum sCoefNames {PERP_DIFFUSION_COEF = 0, PARA_DIFFUSION_COEF,
                    SOURCE_COEF, NUM_SCALAR_COEFS
                   };

   ElectronTotalEnergyCoefs();
};

class CommonCoefs : public EqnCoefficients
{
public:
   enum sCoefNames {IONIZATION_COEF = 0, RECOMBINATION_COEF,
                    CHARGE_EXCHANGE_COEF,
                    NUM_SCALAR_COEFS
                   };
   enum vCoefNames {MAGNETIC_FIELD_COEF = 0, NUM_VECTOR_COEFS};

   CommonCoefs();
};

typedef NeutralDensityCoefs NDCoefs;
typedef IonDensityCoefs IDCoefs;
typedef IonMomentumCoefs IMCoefs;
typedef IonStaticPressureCoefs ISPCoefs;
typedef ElectronStaticPressureCoefs ESPCoefs;
typedef IonTotalEnergyCoefs ITECoefs;
typedef ElectronTotalEnergyCoefs ETECoefs;
typedef CommonCoefs CmnCoefs;

class TransportCoefs
{
private:
   int neqn_;
   Array<EqnCoefficients *> eqnCoefs_;

   void ReadCoefs(common::CoefFactory &cf, std::istream &input);

public:
   TransportCoefs(int neqn)
      : neqn_(neqn),
        eqnCoefs_(neqn+1)
   {
      eqnCoefs_ = NULL;
      eqnCoefs_[0] = new NeutralDensityCoefs;
      eqnCoefs_[1] = new IonDensityCoefs;
      eqnCoefs_[2] = new IonMomentumCoefs;
      // eqnCoefs_[3] = new IonStaticPressureCoefs;
      // eqnCoefs_[4] = new ElectronStaticPressureCoefs;
      eqnCoefs_[3] = new IonTotalEnergyCoefs;
      eqnCoefs_[4] = new ElectronTotalEnergyCoefs;
      eqnCoefs_[5] = new CommonCoefs;
   }

   TransportCoefs(int neqn, common::CoefFactory &cf, std::istream &input);

   ~TransportCoefs()
   {
      for (int i=0; i<=neqn_; i++)
      {
         delete eqnCoefs_[i];
      }
   }

   void LoadCoefs(common::CoefFactory &cf, std::istream &input)
   { ReadCoefs(cf, input); }

   EqnCoefficients & operator()(int i) { return *eqnCoefs_[i]; }
   const EqnCoefficients & operator()(int i) const { return *eqnCoefs_[i]; }

   EqnCoefficients & operator[](int i) { return *eqnCoefs_[i]; }
   const EqnCoefficients & operator[](int i) const { return *eqnCoefs_[i]; }

   NDCoefs & GetNeutralDensityCoefs()
   { return dynamic_cast<NDCoefs&>(*eqnCoefs_[0]); }
   const NDCoefs & GetNeutralDensityCoefs() const
   { return dynamic_cast<const NDCoefs&>(*eqnCoefs_[0]); }

   IDCoefs & GetIonDensityCoefs()
   { return dynamic_cast<IDCoefs&>(*eqnCoefs_[1]); }
   const IDCoefs & GetIonDensityCoefs() const
   { return dynamic_cast<const IDCoefs&>(*eqnCoefs_[1]); }

   IMCoefs & GetIonMomentumCoefs()
   { return dynamic_cast<IMCoefs&>(*eqnCoefs_[2]); }
   const IMCoefs & GetIonMomentumCoefs() const
   { return dynamic_cast<const IMCoefs&>(*eqnCoefs_[2]); }
   /*
   ISPCoefs & GetIonStaticPressureCoefs()
   { return dynamic_cast<ISPCoefs&>(*eqnCoefs_[3]); }
   const ISPCoefs & GetIonStaticPressureCoefs() const
   { return dynamic_cast<const ISPCoefs&>(*eqnCoefs_[3]); }

   ESPCoefs & GetElectronStaticPressureCoefs()
   { return dynamic_cast<ESPCoefs&>(*eqnCoefs_[4]); }
   const ESPCoefs & GetElectronStaticPressureCoefs() const
   { return dynamic_cast<const ESPCoefs&>(*eqnCoefs_[4]); }
   */
   ITECoefs & GetIonTotalEnergyCoefs()
   { return dynamic_cast<ITECoefs&>(*eqnCoefs_[3]); }
   const ITECoefs & GetIonTotalEnergyCoefs() const
   { return dynamic_cast<const ITECoefs&>(*eqnCoefs_[3]); }

   ETECoefs & GetElectronTotalEnergyCoefs()
   { return dynamic_cast<ETECoefs&>(*eqnCoefs_[4]); }
   const ETECoefs & GetElectronTotalEnergyCoefs() const
   { return dynamic_cast<const ETECoefs&>(*eqnCoefs_[4]); }

   CommonCoefs & GetCommonCoefs()
   { return dynamic_cast<CommonCoefs&>(*eqnCoefs_[5]); }
   const CommonCoefs & GetCommonCoefs() const
   { return dynamic_cast<const CommonCoefs&>(*eqnCoefs_[5]); }
};

/*
class ElementSkewGridFunction : public ParGridFunction
{
private:
   ParMesh &pmesh;
   common::L2_FESpace fes;

   inline double cot(double th)
   {
      double s = sin(th);
      if (fabs(s) < 1e-14) { return 1e14; }
      return cos(th) / s;
   }

   void computeSkew();

public:
   ElementSkewGridFunction(ParMesh &_pmesh);

   virtual void Update();
};
*/
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


struct DGParams
{
   double sigma;
   double kappa;
};

struct PlasmaParams
{
   double m_n_amu;
   double m_n_kg;
   double v_n_avg_m_per_s; // Magnitude of average neutral velocity
   double v_n_bar_m_per_s; // Average neutral speed
   double T_n_eV;
   double m_i_amu;
   double m_i_kg;
   int    z_i;
};

class DGAdvectionDiffusionTDO : public TimeDependentOperator
{
private:
   const DGParams & dg_;

   bool imex_;
   int logging_;
   std::string log_prefix_;
   double dt_;

   ParFiniteElementSpace *fes_;
   ParGridFunctionArray  *pgf_;

   Coefficient       *CCoef_;    // Scalar coefficient in front of du/dt
   VectorCoefficient *VCoef_;    // Velocity coefficient
   Coefficient       *dCoef_;    // Scalar diffusion coefficient
   MatrixCoefficient *DCoef_;    // Tensor diffusion coefficient
   Coefficient       *SCoef_;    // Source coefficient

   ScalarVectorProductCoefficient *negVCoef_;   // -1  * VCoef
   ScalarVectorProductCoefficient *dtNegVCoef_; // -dt * VCoef
   ProductCoefficient             *dtdCoef_;    //  dt * dCoef
   ScalarMatrixProductCoefficient *dtDCoef_;    //  dt * DCoef

   Array<int>   dbcAttr_;
   Coefficient *dbcCoef_; // Dirichlet BC coefficient

   Array<int>   nbcAttr_;
   Coefficient *nbcCoef_; // Neumann BC coefficient

   ParBilinearForm  m_;
   ParBilinearForm *a_;
   ParBilinearForm *b_;
   ParBilinearForm *s_;
   ParBilinearForm *k_;
   ParLinearForm   *q_exp_;
   ParLinearForm   *q_imp_;

   HypreParMatrix * M_;
   HypreSmoother M_prec_;
   CGSolver M_solver_;

   // HypreParMatrix * B_;
   // HypreParMatrix * S_;

   mutable ParLinearForm rhs_;
   mutable Vector RHS_;
   mutable Vector X_;

   void initM();
   void initA();
   void initB();
   void initS();
   void initK();
   void initQ();

public:
   DGAdvectionDiffusionTDO(const DGParams & dg,
                           ParFiniteElementSpace &fes,
                           ParGridFunctionArray &pgf,
                           Coefficient &CCoef, bool imex = true);

   ~DGAdvectionDiffusionTDO();

   void SetTime(const double _t);

   void SetLogging(int logging, const std::string & prefix = "");

   void SetAdvectionCoefficient(VectorCoefficient &VCoef);
   void SetDiffusionCoefficient(Coefficient &dCoef);
   void SetDiffusionCoefficient(MatrixCoefficient &DCoef);
   void SetSourceCoefficient(Coefficient &SCoef);

   void SetDirichletBC(Array<int> &dbc_attr, Coefficient &dbc);
   void SetNeumannBC(Array<int> &nbc_attr, Coefficient &nbc);

   virtual void ExplicitMult(const Vector &x, Vector &y) const;
   virtual void ImplicitSolve(const double dt, const Vector &u, Vector &dudt);

   void Update();
};

struct TransPrecParams
{
   int type;
   int log_lvl;

#ifdef MFEM_USE_SUPERLU
   bool l_use_superlu = false;
#endif
#ifdef MFEM_USE_MUMPS
   bool l_use_mumps = false;
#endif
#ifdef MFEM_USE_STRUMPACK
   bool l_use_strumpack = false;
#endif
   bool l_use_algebraic_D_cg = false;
   bool l_use_lor_cg = true;
   bool l_use_air_cg = true;
   bool l_use_schwarz = false;

   bool r_diag_prec = false;
};

struct SolverParams
{
   // Linear solver tolerances
   double lin_abs_tol;
   double lin_rel_tol;
   int lin_max_iter;
   int lin_log_lvl;

   // Newton Solver tolerances
   double newt_abs_tol;
   double newt_rel_tol;
   int newt_max_iter;
   int newt_log_lvl;

   // Steady State tolerances
   double ss_abs_tol;
   double ss_rel_tol;

   TransPrecParams prec;
};

struct CG2DG : Operator
{
   const ParFiniteElementSpace &fes_dg;
   H1_FECollection fec_cg;
   ParFiniteElementSpace fes_cg;
   const HypreParMatrix *P;
   SparseMatrix C;
   mutable Vector z;
   CG2DG(const ParFiniteElementSpace &fes_dg, const Array<int> &cg_ess_tdof_list);
   void Mult(const Vector &x, Vector &y) const;
   void MultTranspose(const Vector &x, Vector &y) const;
   HypreParMatrix *ParallelAssemble(); // Caller must delete returned matrix
};

struct DiscontPSCPreconditioner : Solver
{
   const CG2DG &cg2dg;
   const Solver &cg_solver;
   const Solver &smoother;

   mutable Vector x_z, b_cg, x_cg;
   mutable Vector x_sm;

   DiscontPSCPreconditioner(const CG2DG &cg2dg_,
                            const Solver &cg_solver_,
                            const Solver &smoother_);
   virtual void Mult(const Vector &b, Vector &x) const;
   virtual void SetOperator(const Operator &op);
};

struct AdditivePreconditioner : Solver
{
   const Operator *A;
   const Solver &P1;
   const Solver &P2;

   mutable Vector v;

   AdditivePreconditioner(const Solver &P1_,
                          const Solver &P2_);

   virtual void Mult(const Vector &b, Vector &x) const;
   virtual void SetOperator(const Operator &op);
};

struct MultiplicativePreconditioner : Solver
{
   const Operator *A;
   const Solver &P1;
   const Solver &P2;

   mutable Vector r, v;

   MultiplicativePreconditioner(const Solver &P1_,
                                const Solver &P2_);

   virtual void Mult(const Vector &b, Vector &x) const;
   virtual void SetOperator(const Operator &op);
};

/** The DGTransportTDO class is designed to be used with an implicit
    ODESolver to solve a specific set of coupled transport equations
    using a Discontinuous Galerkin (DG) discretization of the relevant
    PDEs.

    This system of transport equations consists of mass conservation
    equations for one species of neutrals and one species of ions,
    momentum conservation for the ion species, and energy equations
    for ions and electrons.

    The field variables are density of neutrals, density of ions,
    velocity of ions in a direction parallel to the magnetic field,
    and the ion and electron temperatures. The electron density is
    assumed to depend on the ion density in such a way that
    quasi-neutrality is maintained e.g. n_e = z_i n_i.

    The ODE solvers in MFEM integrate, in time, equations of the form
       du/dt = f(u,t)

    Where f(u,t) can be a non-linear function of the field quantities,
    u(x,t). The primary function of this class is to implement the
    ImplicitSolve method which finds a du/dt which satisfies
       du/dt = f(u + dt * du/dt, t + dt)
    Since f(u,t) is generally non-linear this is computed using a
    Newton solver.
*/
class DGTransportTDO : public TimeDependentOperator
{
   // friend class TransportPrec;
private:
   int logging_;
   int op_flag_;

   ParFiniteElementSpace &fes_;
   ParFiniteElementSpace &ffes_;
   ParGridFunctionArray  &yGF_;
   ParGridFunctionArray  &kGF_;

   Array<int> &offsets_;

   SolverParams tol_;

   Vector kMax_;
   Array<bool> ss_;

   // Data collection used to write data files
   DataCollection * dc_;

   // Sockets used to communicate with GLVis
   std::map<std::string, socketstream*> socks_;

   class NLOperator : public Operator
   {
   protected:
      const DGParams &dg_;

      int logging_;
      std::string log_prefix_;

      int index_;
      std::string eqn_name_;
      std::string field_name_;
      double dt_;
      ParFiniteElementSpace &fes_;
      ParMesh               &pmesh_;
      ParGridFunctionArray  &yGF_;
      ParGridFunctionArray  &kGF_;

      VectorCoefficient & B3Coef_;

      Array<StateVariableGridFunctionCoef*>  yCoefPtrs_;
      Array<StateVariableGridFunctionCoef*>  kCoefPtrs_;
      Array<StateVariableSumCoef*>          ykCoefPtrs_;

      mutable Array<int> vdofs_;
      mutable Array<int> vdofs2_;
      mutable DenseMatrix elmat_;
      mutable DenseMatrix elmat_k_;
      mutable Vector elvec_;
      mutable Vector locvec_;
      mutable Vector locdvec_;

      typedef std::vector<Array<BilinearFormIntegrator*> > vArrayBFI;

      // Domain integrators for time derivatives of field variables
      vArrayBFI dbfi_m_;  // Domain Integrators
      // Array<Array<StateVariableCoef*> >      dbfi_mc_; // Domain Integrators

      // Domain integrators for field variables at next time step
      vArrayBFI dbfi_;  // Domain Integrators
      vArrayBFI fbfi_;  // Interior Face Integrators
      vArrayBFI bfbfi_; // Boundary Face Integrators
      std::vector<Array<Array<int>*> > bfbfi_marker_; ///< Entries are owned.

      // Domain integrators for source terms
      Array<LinearFormIntegrator*> dlfi_;  // Domain Integrators
      Array<LinearFormIntegrator*> bflfi_; // Boundary Face Integrators
      Array<Array<int>*>           bflfi_marker_; ///< Entries are owned.

      Array<ParBilinearForm*> blf_; // Bilinear Form Objects for Gradients
      Array<ParBilinearForm*> cgblf_; // Bilinear Form Objects for Gradients

      HypreParMatrix * Dmat_ = NULL;
      CG2DG *cg2dg_ = NULL;
      HypreParMatrix *CG2DGmat_ = NULL;
      Solver *D_amg_ = NULL;
      ParLORDiscretization *D_lor_ = NULL;
      HypreSmoother *D_smoother_ = NULL;
      HypreParMatrix *D_cg_ = NULL;
      Solver *dg_precond_ = NULL;

      Solver *D_mult_ = NULL;
      Solver *D_schwarz_ = NULL;
      Vector D_diag_;

      Array<int> cg_ess_tdof_list;
      /*
      bool use_algebraic_D_cg = false;
      bool use_lor_cg = true;
      bool use_air_cg = true;
      bool use_schwarz = false;
      */
      int term_flag_;
      int vis_flag_;

      // Data collection used to write data files
      DataCollection * dc_;

      // Sockets used to communicate with GLVis
      std::map<std::string, socketstream*> socks_;

      NLOperator(const DGParams & dg,
                 int index,
                 const std::string &eqn_name,
                 const std::string &field_name,
                 ParGridFunctionArray & yGF,
                 ParGridFunctionArray & kGF,
                 VectorCoefficient & B3Coef,
                 int term_flag, int vis_flag, int logging,
                 const std::string & log_prefix);


   public:

      virtual ~NLOperator();

      void SetLogging(int logging, const std::string & prefix = "");

      virtual void SetTimeStep(double dt);

      virtual void Mult(const Vector &k, Vector &y) const;

      virtual void Update();

      virtual void PrepareGradient() { ; }
      virtual Operator *GetGradientBlock(int i);

      virtual Solver *GetPreconditioner(const TransPrecParams &pparams);

      inline bool CheckTermFlag(int flag) { return (term_flag_>> flag) & 1; }

      inline bool CheckVisFlag(int flag) { return (vis_flag_>> flag) & 1; }

      virtual int GetDefaultVisFlag() = 0;

      virtual void RegisterDataFields(DataCollection & dc);

      virtual void PrepareDataFields();

      virtual void InitializeGLVis() = 0;

      virtual void DisplayToGLVis() = 0;
   };

   class TransportOp : public NLOperator
   {
   protected:
      Array<StateVariableCoef*>    svscoefs_;
      Array<StateVariableVecCoef*> svvcoefs_;
      Array<StateVariableMatCoef*> svmcoefs_;
      Array<ProductCoefficient*>             dtSCoefs_;
      Array<ProductCoefficient*>             negdtSCoefs_;
      Array<ScalarVectorProductCoefficient*> dtVCoefs_;
      Array<ScalarVectorProductCoefficient*> negdtVCoefs_;
      Array<ScalarMatrixProductCoefficient*> dtMCoefs_;
      Array<Coefficient*>       sCoefs_;
      Array<VectorCoefficient*> vCoefs_;
      Array<MatrixCoefficient*> mCoefs_;
      std::vector<socketstream*> sout_;
      ParGridFunction coefGF_;

      ParFiniteElementSpace *vfes_;
      ParFiniteElementSpace *h1_fes_;

      const PlasmaParams &plasma_;

      double m_n_kg_;
      double T_n_eV_;
      double v_n_;
      double m_i_kg_;
      int    z_i_;
      /*
      StateVariableGridFunctionCoef &nn0Coef_;
      StateVariableGridFunctionCoef &ni0Coef_;
      StateVariableGridFunctionCoef &vi0Coef_;
      StateVariableGridFunctionCoef &Ti0Coef_;
      StateVariableGridFunctionCoef &Te0Coef_;
      */
      StateVariableCoef &nnCoef_;
      StateVariableCoef &niCoef_;
      StateVariableCoef &viCoef_;
      StateVariableCoef &TiCoef_;
      StateVariableCoef &TeCoef_;

      StateVariableConstantCoef vnAvgCoef_; // Magnitude of average velocity
      StateVariableConstantCoef vnBarCoef_; // Average speed of neutrals
      StateVariableConstantCoef TnCoef_;
      StateVariableConstantCoef ziCoef_;
      StateVariableProductCoef  neCoef_;

      StateVariableGridFunctionCoef &dTe0Coef_;

      CoulombLogEICoef lnLambda_;

      const AdvectionDiffusionBC & bcs_;
      const CoupledBCs & cbcs_;

      const CommonCoefs & cmncoefs_;

      GridFunctionCoefficient elOrdCoef_;
      GridFunctionCoefficient hCoef_;

      VectorCoefficient & B3Coef_;
      VectorXYCoefficient BxyCoef_;
      UnitVectorXYCoefficient UnitBxyCoef_;

      Coefficient       * massCoef_;
      Coefficient       * diffusionCoef_;
      MatrixCoefficient * diffusionMatrixCoef_;

      ApproxIonizationRate    izCoef_;
      ApproxRecombinationRate rcCoef_;
      ApproxChargeExchangeRate cxCoef_;

      IonizationSourceCoef    SizDefCoef_;
      RecombinationSinkCoef   SrcDefCoef_;
      ChargeExchangeSinkCoef  ScxDefCoef_;

      StateVariableCoef & SizCoef_;
      StateVariableCoef & SrcCoef_;
      StateVariableCoef & ScxCoef_;

      ParFluxVectors flux_vis_;

      TransportOp(const DGParams & dg,
                  const PlasmaParams & plasma, int index,
                  const std::string &eqn_name,
                  const std::string &field_name,
                  ParFiniteElementSpace * vfes,
                  ParFiniteElementSpace * h1_fes,
                  ParGridFunctionArray & yGF,
                  ParGridFunctionArray & kGF,
                  ParGridFunction & elOrdGF,
                  ParGridFunction & hGF,
                  const AdvectionDiffusionBC & bcs,
                  const CoupledBCs & cbcs,
                  const CommonCoefs & common_coefs,
                  VectorCoefficient & B3Coef,
                  int term_flag, int vis_flag,
                  int logging,
                  const std::string & log_prefix);

      /** Sets the time derivative on the left hand side of the equation to be:
             d MCoef / dt
      */
      void SetTimeDerivativeTerm(StateVariableCoef &MCoef);

      /** Sets the diffusion term on the right hand side of the equation
          to be:
             Div(DCoef Grad y[index])
          where index is the index of the equation.
       */
      void SetDiffusionTerm(StateVariableCoef &DCoef);
      void SetDiffusionTerm(StateVariableMatCoef &DCoef);
      /*
       void SetAnisoDiffusionTerm(StateVariableMatCoef &DCoef,
                                  Coefficient &SkewCoef,
                                  double D_min, double D_max);
      */
      /** Sets the advection-diffusion term on the right hand side of the
          equation to be:
             Div(DCoef Grad y[index] - VCoef y[index])
           where index is the index of the equation.
       */
      void SetAnisotropicDiffusionTerm(StateVariableMatCoef &DCoef,
                                       Coefficient *DParaCoef,
                                       Coefficient *DPerpCoef);

      /** Sets the advection-diffusion term on the right hand side of the
          equation to be:
             Div(DCoef Grad y[index] - VCoef y[index])
           where index is the index of the equation.
       */
      // void SetAdvectionDiffusionTerm(StateVariableMatCoef &DCoef,
      //                              StateVariableVecCoef &VCoef,
      //                              Coefficient *DParaCoef,
      //                              Coefficient *DPerpCoef);
      void SetDiffusionTermGradient(StateVariableMatCoef &DCoef);

      /** Sets the advection term on the right hand side of the
      equation to be:
             Div(VCoef y[index])
          where index is the index of the equation.
       */
      // void SetAdvectionTerm(StateVariableVecCoef &VCoef/*, bool bc = false*/);

      void SetSourceTerm(StateVariableCoef &SCoef, double s = 1.0);
      void SetSourceTermGradient(StateVariableCoef &SCoef, double s = 1.0);
      void SetBdrSourceTerm(StateVariableCoef &SCoef,
                            StateVariableVecCoef &VCoef);

      void SetOutflowBdrTerm(StateVariableVecCoef &VCoef,
                             const Array<CoefficientByAttr*> & obc);
      void SetRecyclingBdrSourceTerm(const RecyclingBC & rbc);

   public:
      virtual ~TransportOp();

      virtual void Update();

      virtual void SetTime(double t);
      virtual void SetTimeStep(double dt);

      virtual void InitializeGLVis();
      virtual void DisplayToGLVis();

      inline Coefficient       * GetMassCoef() { return massCoef_; }
      inline Coefficient       * GetDiffusionCoef() { return diffusionCoef_; }
      inline MatrixCoefficient * GetDiffusionMatrixCoef()
      { return diffusionMatrixCoef_; }
   };

   class AdvTransportOp : public TransportOp
   {
   protected:

      VectorCoefficient     * advectionCoef_;
      common::L2_ParFESpace * l2_fes_0_;
      common::H1_ParFESpace * h1_fes_1_;
      ParGridFunction       * elOrdDiscGF_;
      ParGridFunction       * elOrdContGF_;
      ParGridFunction       * OscDiscGF_;
      ParGridFunction       * OscContGF_;
      ParGridFunction       * hDiscGF_;
      ParGridFunction       * hContGF_;
      GridFunctionCoefficient elOrdDiscCoef_;
      GridFunctionCoefficient elOrdContCoef_;
      GridFunctionCoefficient OscDiscCoef_;
      GridFunctionCoefficient OscContCoef_;
      GridFunctionCoefficient hDiscCoef_;
      GridFunctionCoefficient hContCoef_;
      SoundSpeedCoef          CsCoef_;

      AdvTransportOp(const DGParams & dg,
                     const PlasmaParams & plasma, int index,
                     const std::string &eqn_name,
                     const std::string &field_name,
                     ParFiniteElementSpace * vfes,
                     ParFiniteElementSpace * h1_fes,
                     ParGridFunctionArray & yGF,
                     ParGridFunctionArray & kGF,
                     ParGridFunction & elOrdGF,
                     ParGridFunction & hGF,
                     const AdvectionDiffusionBC & bcs,
                     const CoupledBCs & cbcs,
                     const CommonCoefs & common_coefs,
                     VectorCoefficient & B3Coef,
                     int term_flag, int vis_flag,
                     int logging,
                     const std::string & log_prefix);

      /** Sets the advection-diffusion term on the right hand side of the
          equation to be:
             Div(DCoef Grad y[index] - VCoef y[index])
           where index is the index of the equation.
       */
      void SetAdvectionDiffusionTerm(StateVariableMatCoef &DCoef,
                                     StateVariableVecCoef &VCoef,
                                     Coefficient *DParaCoef,
                                     Coefficient *DPerpCoef);

      /** Sets the advection term on the right hand side of the
      equation to be:
             Div(VCoef y[index])
          where index is the index of the equation.
       */
      void SetAdvectionTerm(StateVariableVecCoef &VCoef/*, bool bc = false*/);

   public:
      virtual ~AdvTransportOp();

      virtual void SetTime(double t);
      // virtual void SetTimeStep(double dt);

      virtual void Update();

      // virtual void InitializeGLVis();
      // virtual void DisplayToGLVis();

      inline VectorCoefficient * GetAdvectionCoef() { return advectionCoef_; }
   };

   /** The NeutralDensityOp is an mfem::Operator designed to work with
       a NewtonSolver as one row in a block system of non-linear
       transport equations.  Specifically, this operator models the
       mass conservation equation for a neutral species.

          d n_n / dt = Div(D_n Grad(n_n)) + S_n

       Where the diffusion coefficient D_n is a function of n_e and T_e
       (the electron density and temperature respectively) and the
       source term S_n is a function of n_e, T_e, and n_n.  Note that n_e is
       not a state variable but is related to n_i by the simple relation
          n_e = z_i n_i
       where z_i is the charge of the ions and n_i is the ion density.

       To advance this equation in time we need to find k_nn = d n_n / dt
       which satisfies:
          k_nn - Div(D_n(n_e, T_e) Grad(n_n + dt k_nn))
               - S_n(n_e, T_e, n_n + dt k_nn) = 0
       Where n_e and T_e are also evaluated at the next time step.  This is
       done with a Newton solver which needs the Jacobian of this block of
       equations.

       The diagonal block is given by:
          1 - dt Div(D_n Grad) - dt d S_n / d n_n

       The other non-trivial blocks are:
          - dt Div(d D_n / d n_i Grad(n_n)) - dt d S_n / d n_i
          - dt Div(d D_n / d T_e Grad(n_n)) - dt d S_n / d T_e

       The blocks of the Jacobian will be assembled finite element matrices.
       For the diagonal block we need a mass integrator with coefficient
       (1 - dt d S_n / d n_n), and a set of integrators to model the DG
       diffusion operator with coefficient (dt D_n).

       The off-diagonal blocks will consist of a mass integrator with
       coefficient (-dt d S_n / d n_i) or (-dt d S_n / d
       T_e). Currently, (-dt d S_n / d T_e) is not implemented.
    */
   class NeutralDensityOp : public TransportOp
   {
   private:
      enum TermFlag {DIFFUSION_TERM = 0,
                     RECOMBINATION_SOURCE_TERM,
                     IONIZATION_SINK_TERM,
                     SOURCE_TERM,
                     RECYCLING_BDR_SOURCE_TERM,
                     NUM_TERMS
                    };
      enum VisField {DIFFUSION_COEF = 1,
                     RECOMBINATION_SOURCE_COEF,
                     IONIZATION_SINK_COEF,
                     SOURCE_COEF,
                     DIFFUSIVE_FLUX,
                     NUM_FIELDS_PLUS_ONE
                    };

      const NeutralDensityCoefs & ndcoefs_;

      NeutralDiffusionCoef      DDefCoef_; // Default diffusion coef
      StateVariableStandardCoef DCoef_;

      ParGridFunction * DGF_;
      ParGridFunction * SrcGF_;
      ParGridFunction * SizGF_;
      ParGridFunction * SGF_;

   public:
      friend class NeutralDensityTerms;

      NeutralDensityOp(const DGParams & dg,
                       const PlasmaParams & plasma,
                       ParFiniteElementSpace & h1_fes,
                       ParGridFunctionArray & yGF,
                       ParGridFunctionArray & kGF,
                       ParGridFunction & elOrdGF,
                       ParGridFunction & hGF,
                       const AdvectionDiffusionBC & bcs,
                       const CoupledBCs & cbcs,
                       const NDCoefs & ndcoefs,
                       const CommonCoefs & cmncoefs,
                       VectorCoefficient & B3Coef,
                       int term_flag,
                       int vis_flag, int logging,
                       const std::string & log_prefix);

      virtual ~NeutralDensityOp();

      virtual void SetTime(double t);
      virtual void SetTimeStep(double dt);

      void Update();

      virtual int GetDefaultVisFlag() { return 7; }

      virtual void RegisterDataFields(DataCollection & dc);

      virtual void PrepareDataFields();
   };

   /** The IonDensityOp is an mfem::Operator designed to worth with a
       NewtonSolver as one row in a block system of non-linear
       transport equations.  Specifically, this operator models the
       mass conservation equation for a single ion species.

       d n_i / dt = Div(D_i Grad n_i)) - Div(v_i n_i b_hat) + S_i

       Where the diffusion coefficient D_i is a function of the
       magnetic field direction, v_i is the velocity of the ions
       parallel to B, and the source term S_i is a function of the
       electron and neutral densities as well as the electron
       temperature.

       To advance this equation in time we need to find k_ni = d n_i / dt
       which satisfies:
          k_ni - Div(D_i Grad(n_i + dt k_ni)) + Div(v_i (n_i + dt k_ni) b_hat)
               - S_i(n_e + z_i dt k_ni, T_e, n_n) = 0
       Where n_n and T_e are also evaluated at the next time step.  This is
       done with a Newton solver which needs the Jacobian of this block of
       equations.

       The diagonal block is given by:
          1 - dt Div(D_i Grad) + dt Div(v_i b_hat) - dt d S_i / d n_i

       The other non-trivial blocks are:
          - dt d S_i / d n_n
          + dt Div(n_i b_hat)
          - dt d S_i / d T_e

       The blocks of the Jacobian will be assembled finite element
       matrices.  For the diagonal block we need a mass integrator
       with coefficient (1 - dt d S_i / d n_i), a set of integrators
       to model the DG diffusion operator with coefficient (dt D_i),
       and a weak divergence integrator with coefficient (dt v_i).

       The off-diagonal blocks will consist of a mass integrator with
       coefficient (-dt d S_i / d n_n) or (-dt d S_i / d T_e).
       Currently, (dt Div(n_i b_hat)) and (-dt d S_i / d T_e) are not
       implemented.
    */
   class IonDensityOp : public AdvTransportOp
   {
   private:
      enum TermFlag {DIFFUSION_TERM = 0,
                     ADVECTION_TERM,
                     IONIZATION_SOURCE_TERM,
                     RECOMBINATION_SINK_TERM,
                     SOURCE_TERM,
                     RECYCLING_BDR_SINK_TERM,
                     NUM_TERMS
                    };
      enum VisField {DIFFUSION_PARA_COEF = 1,
                     DIFFUSION_PERP_COEF,
                     ADVECTION_COEF,
                     IONIZATION_SOURCE_COEF,
                     RECOMBINATION_SINK_COEF,
                     SOURCE_COEF,
                     DIFFUSIVE_FLUX,
                     ADVECTIVE_FLUX,
                     NUM_FIELDS_PLUS_ONE
                    };

      const IonDensityCoefs & idcoefs_;

      const ArtViscParams & av_;

      StateVariableConstantCoef   DPerpConstCoef_;
      IonDensityParaDiffusionCoef DParaCoef_;
      StateVariableCoef *         DParaCoefPtr_;
      StateVariableCoef *         DPerpCoefPtr_;
      Aniso2DDiffusionCoef        DCoef_;

      IonAdvectionCoef        ViCoef_;

      ParGridFunction * DParaGF_;
      ParGridFunction * DPerpGF_;
      ParGridFunction * AdvGF_;
      ParGridFunction * SizGF_;
      ParGridFunction * SrcGF_;
      ParGridFunction * SGF_;

   public:
      friend class IonDensityTerms;

      IonDensityOp(const DGParams & dg,
                   const ArtViscParams & av,
                   const PlasmaParams & plasma,
                   ParFiniteElementSpace & vfes,
                   ParFiniteElementSpace & h1_fes,
                   ParGridFunctionArray & yGF,
                   ParGridFunctionArray & kGF,
                   ParGridFunction & elOrdGF,
                   ParGridFunction & hGF,
                   const AdvectionDiffusionBC & bcs,
                   const CoupledBCs & cbcs,
                   const IDCoefs & idcoefs,
                   const CommonCoefs & cmncoefs,
                   VectorCoefficient & B3Coef,
                   double DPerp,
                   int term_flag, int vis_flag, int logging,
                   const std::string & log_prefix);

      virtual ~IonDensityOp();

      virtual void SetTime(double t);
      virtual void SetTimeStep(double dt);

      void Update();

      virtual int GetDefaultVisFlag() { return 27; }

      virtual void RegisterDataFields(DataCollection & dc);

      virtual void PrepareDataFields();
   };

   /** The IonMomentumOp is an mfem::Operator designed to work with a
       NewtonSolver as one row in a block system of non-linear
       transport equations.  Specifically, this operator models the
       momentum conservation equation for a single ion species.

       m_i v_i d n_i / dt + m_i n_i d v_i / dt
          = Div(Eta Grad v_i) - Div(m_i n_i v_i v_i) - b.Grad(p_i + p_e)
          + m_i v_n S_iz - m_i v_i S_rc + m_i (v_n - v_i) S_cx

       Where the diffusion coefficient Eta is a function of the
       magnetic field, ion density, and ion temperature.

       To advance this equation in time we need to find k_vi = d v_i / dt
       which satisfies:
          m_i n_i k_vi - Div(Eta Grad(v_i + dt k_vi))
             + Div(m_i n_i v_i (v_i + dt k_vi)) + b.Grad(p_i + p_e)
             - m_i v_n S_iz + m_i (v_i + dt k_vi) S_rc
             - m_i (v_n - v_i - dt k_vi) S_cx = 0
       Where n_i, p_i, and p_e are also evaluated at the next time
       step.  This is done with a Newton solver which needs the
       Jacobian of this block of equations.

       The diagonal block is given by:
          m_i n_i - dt Div(Eta Grad) + dt Div(m_i n_i v_i)
             + dt m_i (S_rc + S_cx)
       MLS: Why is the advection term not doubled?

       The other non-trivial blocks are:
          m_i v_i - dt Div(d Eta / d n_i Grad(v_i)) + dt Div(m_i v_i v_i)
          - dt Div(d Eta / d T_i Grad(v_i)) + dt b.Grad(d p_i / d T_i)
          + dt b.Grad(d p_e / d T_e)

       Currently, the static pressure terms and the derivatives of Eta
       do not contribute to the Jacobian.
    */
   class IonMomentumOp : public AdvTransportOp
   {
   private:
      enum TermFlag {DIFFUSION_TERM = 0,
                     ADVECTION_TERM, GRADP_SOURCE_TERM,
                     DIVBP_SOURCE_TERM,
                     IONIZATION_SOURCE_TERM,
                     RECOMBINATION_SINK_TERM, CHARGE_EXCHANGE_SOURCE_TERM,
                     SOURCE_TERM,
                     NUM_TERMS
                    };
      enum VisField {DIFFUSION_PARA_COEF = 1, DIFFUSION_PERP_COEF,
                     ADVECTION_COEF, GRADP_SOURCE_COEF,
                     IONIZATION_SOURCE_COEF,
                     RECOMBINATION_SINK_COEF, CHARGE_EXCHANGE_SOURCE_COEF,
                     SOURCE_COEF,
                     ION_PARA_MOMENTUM,
                     NUM_FIELDS_PLUS_ONE
                    };

      const IonMomentumCoefs & imcoefs_;

      const ArtViscParams & av_;

      /*
      common::L2_ParFESpace * l2_fes_0_;
      common::H1_ParFESpace * h1_fes_1_;
      ParGridFunction       * elOrdDiscGF_;
      ParGridFunction       * elOrdContGF_;
      ParGridFunction       * OscDiscGF_;
      ParGridFunction       * OscContGF_;
      ParGridFunction       * hDiscGF_;
      ParGridFunction       * hContGF_;
      GridFunctionCoefficient elOrdDiscCoef_;
      GridFunctionCoefficient elOrdContCoef_;
      GridFunctionCoefficient OscDiscCoef_;
      GridFunctionCoefficient OscContCoef_;
      GridFunctionCoefficient hDiscCoef_;
      GridFunctionCoefficient hContCoef_;
      */
      // SoundSpeedCoef          CsCoef_;

      double DPerpConst_;
      StateVariableConstantCoef DPerpCoef_;

      IonMomentumParaCoef            momCoef_;
      IonMomentumParaDiffusionCoef   EtaParaCoef_;
      IonMomentumPerpDiffusionCoef   EtaPerpCoef_;
      StateVariableCoef *            EtaParaCoefPtr_;
      StateVariableCoef *            EtaPerpCoefPtr_;
      Aniso2DDiffusionCoef           EtaCoef_;

      IonMomentumAdvectionCoef miniViCoef_;

      NegGradPressureCoefficient negGradPCoef_;
      NegBPressureCoefficient negBPCoef_;
      NegPressureCoefficient negPCoef_;

      IonMomentumIonizationCoef     SIZCoef_;
      IonMomentumRecombinationCoef  SRCCoef_;
      IonMomentumChargeExchangeCoef SCXCoef_;

      ParGridFunction * EtaParaGF_;
      ParGridFunction * EtaPerpGF_;
      ParGridFunction * AdvGF_;
      ParGridFunction * MomParaGF_;
      ParGridFunction * SGPGF_;
      ParGridFunction * SIZGF_;
      ParGridFunction * SRCGF_;
      ParGridFunction * SCXGF_;
      ParGridFunction * SGF_;

      /** Sets a divergence term on the right hand side of the
      equation to be:
             Div(VCoef)

       */
      void SetDivergenceTerm(StateVariableVecCoef &VCoef);
      void SetDivergenceTerm(StateVariableCoef &Coef,
                             VectorCoefficient &VCoef);

   public:
      friend class IonMomentumTerms;

      IonMomentumOp(const DGParams & dg,
                    const ArtViscParams & av,
                    const PlasmaParams & plasma,
                    ParFiniteElementSpace & vfes,
                    ParFiniteElementSpace & h1_fes,
                    ParGridFunctionArray & yGF, ParGridFunctionArray & kGF,
                    ParGridFunction & elOrdGF,
                    ParGridFunction & hGF,
                    const AdvectionDiffusionBC & bcs,
                    const CoupledBCs & cbcs,
                    const IMCoefs & imcoefs,
                    const CmnCoefs & cmncoefs,
                    VectorCoefficient & B3Coef,
                    double DPerp,
                    int term_flag, int vis_flag, int logging,
                    const std::string & log_prefix);

      virtual ~IonMomentumOp();

      virtual void SetTime(double t);
      virtual void SetTimeStep(double dt);

      virtual void Update();
      virtual void PrepareGradient();

      virtual int GetDefaultVisFlag() { return 127; }

      virtual void RegisterDataFields(DataCollection & dc);

      virtual void PrepareDataFields();
   };

   /** The IonStaticPressureOp is an mfem::Operator designed to work
       with a NewtonSolver as one row in a block system of non-linear
       transport equations.  Specifically, this operator models the
       static pressure equation (related to conservation of energy)
       for a single ion species.

       1.5 T_i d n_i / dt + 1.5 n_i d T_i / dt
          = Div(n_i Chi_i Grad(T_i) - 2.5 n_i v_i T_i)

       Where the diffusion coefficient Chi_i is a function of the
       magnetic field direction, ion density and temperature.

       MLS: Clearly this equation is incomplete.  We stopped at this
       point to focus on implementing a non-linear Robin boundary
       condition.

       To advance this equation in time we need to find
       k_Ti = d T_i / dt which satisfies:
          (3/2)(T_i d n_i / dt + n_i k_Ti) - Div(Chi_i Grad T_i) = 0
       Where n_i is also evaluated at the next time step.  This is
       done with a Newton solver which needs the Jacobian of this
       block of equations.

       The diagonal block is given by:
          1.5 n_i - dt Div(Chi_i Grad)

       The other non-trivial blocks are:
          1.5 T_i - dt Div(d Chi_i / d n_i Grad(T_i))

       MLS: Many more terms will arise once the full equation is implemented.
   */
   class IonStaticPressureOp : public AdvTransportOp
   {
   private:
      enum TermFlag {DIFFUSION_TERM = 0, ADVECTION_TERM, SOURCE_TERM, NUM_TERMS};
      enum VisField {DIFFUSION_PARA_COEF = 1, DIFFUSION_PERP_COEF,
                     SOURCE_COEF, NUM_FIELDS_PLUS_ONE
                    };

      const IonStaticPressureCoefs & ispcoefs_;

      double ChiPerpConst_;

      StaticPressureCoef               presCoef_;
      StaticPressureAdvectionCoef      aniViCoef_;
      IonThermalParaDiffusionCoef      ChiParaCoef_;
      StateVariableConstantCoef        ChiPerpCoef_;
      StateVariableCoef *              ChiParaCoefPtr_;
      StateVariableCoef *              ChiPerpCoefPtr_;
      Aniso2DDiffusionCoef             ChiCoef_;
      ProductCoefficient               nChiParaCoef_;
      ProductCoefficient               nChiPerpCoef_;
      StateVariableScalarMatrixProductCoef nChiCoef_;

      ParGridFunction * ChiParaGF_;
      ParGridFunction * ChiPerpGF_;
      ParGridFunction * SGF_;

   public:
      IonStaticPressureOp(const DGParams & dg,
                          const PlasmaParams & plasma,
                          ParGridFunctionArray & yGF,
                          ParGridFunctionArray & kGF,
                          ParGridFunction & elOrdGF,
                          ParGridFunction & hGF,
                          const AdvectionDiffusionBC & bcs,
                          const CoupledBCs & cbcs,
                          const ISPCoefs & ispcoefs,
                          const CmnCoefs & cmncoefs,
                          VectorCoefficient & B3Coef,
                          double ChiPerp,
                          int term_flag, int vis_flag, int logging,
                          const std::string & log_prefix);

      virtual ~IonStaticPressureOp();

      virtual void SetTimeStep(double dt);

      void Update();

      virtual int GetDefaultVisFlag() { return 4; }

      virtual void RegisterDataFields(DataCollection & dc);

      virtual void PrepareDataFields();
   };

   /** The ElectronStaticPressureOp is an mfem::Operator designed to
       work with a NewtonSolver as one row in a block system of
       non-linear transport equations.  Specifically, this operator
       models the static pressure equation (related to conservation of
       energy) for the flow of electrons.

       1.5 T_e d n_e / dt + 1.5 n_e d T_e / dt
          = Div(Chi_e Grad(T_e) - 2.5 n_e v_i T_e)

       Where the diffusion coefficient Chi_e is a function of the
       magnetic field direction, electron density and temperature.

       MLS: Clearly this equation is incomplete.  We stopped at this
       point to focus on implementing a non-linear Robin boundary
       condition.

       To advance this equation in time we need to find
       k_Te = d T_e / dt which satisfies:
          (3/2)(T_e d n_e / dt + n_e k_Te) - Div(Chi_e Grad T_e) = 0
       Where n_e is also evaluated at the next time step.  This is
       done with a Newton solver which needs the Jacobian of this
       block of equations.

       The diagonal block is given by:
          1.5 n_e - dt Div(Chi_e Grad)

       The other non-trivial blocks are:
          1.5 T_e - dt Div(d Chi_e / d n_e Grad(T_e))

       MLS: Many more terms will arise once the full equation is implemented.
   */
   class ElectronStaticPressureOp : public AdvTransportOp
   {
   private:
      enum TermFlag {DIFFUSION_TERM = 0, ADVECTION_TERM, SOURCE_TERM, NUM_TERMS};
      enum VisField {DIFFUSION_PARA_COEF = 1, DIFFUSION_PERP_COEF,
                     SOURCE_COEF, NUM_FIELDS_PLUS_ONE
                    };

      const ElectronStaticPressureCoefs & espcoefs_;

      double ChiPerpConst_;

      StaticPressureCoef               presCoef_;
      StaticPressureAdvectionCoef      aneViCoef_;
      ElectronThermalParaDiffusionCoef ChiParaCoef_;
      StateVariableConstantCoef        ChiPerpCoef_;
      StateVariableCoef *              ChiParaCoefPtr_;
      StateVariableCoef *              ChiPerpCoefPtr_;
      Aniso2DDiffusionCoef             ChiCoef_;
      ProductCoefficient               nChiParaCoef_;
      ProductCoefficient               nChiPerpCoef_;
      StateVariableScalarMatrixProductCoef nChiCoef_;

      ParGridFunction * ChiParaGF_;
      ParGridFunction * ChiPerpGF_;
      ParGridFunction * SGF_;

   public:
      ElectronStaticPressureOp(const DGParams & dg,
                               const PlasmaParams & plasma,
                               ParGridFunctionArray & yGF,
                               ParGridFunctionArray & kGF,
                               ParGridFunction & elOrdGF,
                               ParGridFunction & hGF,
                               const AdvectionDiffusionBC & bcs,
                               const CoupledBCs & cbcs,
                               const ESPCoefs & espcoefs,
                               const CmnCoefs & cmncoefs,
                               VectorCoefficient & B3Coef,
                               double ChiPerp,
                               int term_flag, int vis_flag,
                               int logging,
                               const std::string & log_prefix);

      virtual ~ElectronStaticPressureOp();

      virtual void SetTimeStep(double dt);

      void RegisterDataFields(DataCollection & dc);

      void PrepareDataFields();

      void Update();

      virtual int GetDefaultVisFlag() { return 4; }
   };

   class TotalEnergyOp : public AdvTransportOp
   {
   protected:

      IonElectronHeatExchangeCoef QiCoef_;

      StateVariableConstantCoef kBCoef_;
      StateVariableConstantCoef phiIZCoef_;
      StateVariableProductCoef kBphiIZCoef_;
      StateVariableStandardVecCoef BSVCoef_;

      TotalEnergyOp(const DGParams & dg,
                    const PlasmaParams & plasma, int index,
                    const std::string &eqn_name,
                    const std::string &field_name,
                    ParFiniteElementSpace & vfes,
                    ParFiniteElementSpace & h1_fes,
                    ParGridFunctionArray & yGF,
                    ParGridFunctionArray & kGF,
                    ParGridFunction & elOrdGF,
                    ParGridFunction & hGF,
                    const AdvectionDiffusionBC & bcs,
                    const CoupledBCs & cbcs,
                    const CommonCoefs & common_coefs,
                    VectorCoefficient & B3Coef,
                    int term_flag, int vis_flag,
                    int logging,
                    const std::string & log_prefix);

   public:

      void SetKineticEnergyAdvectionTerm(StateVariableVecCoef &VCoef);

      virtual void SetTime(double t);
   };

   class IonTotalEnergyOp : public TotalEnergyOp
   {
   private:
      enum TermFlag {DIFFUSION_TERM = 0,
                     ADVECTION_TERM, KE_ADVECTION_TERM,
                     IONIZATION_SOURCE_TERM, RECOMBINATION_SINK_TERM,
                     CHARGE_EXCHANGE_SOURCE_TERM,
                     EQUIPARTITION_SOURCE_TERM, SOURCE_TERM, NUM_TERMS
                    };

      enum VisField {DIFFUSION_PARA_COEF = 1, DIFFUSION_PERP_COEF,
                     ADVECTION_COEF, IONIZATION_SOURCE_COEF,
                     RECOMBINATION_SINK_COEF, CHARGE_EXCHANGE_SOURCE_COEF,
                     EQUIPARTITION_SOURCE_COEF, SOURCE_COEF,
                     ION_TOTAL_ENERGY,
                     DIFFUSIVE_FLUX,
                     ADVECTIVE_FLUX,
                     NUM_FIELDS_PLUS_ONE
                    };

      const IonTotalEnergyCoefs & itecoefs_;

      double ChiPerpConst_;

      TotalEnergyCoef                  totEnergyCoef_;
      KineticEnergyCoef                kinEnergyCoef_;
      TotalEnergyAdvectionCoef         advFluxCoef_;
      StaticPressureAdvectionCoef      aniViCoef_;
      IonThermalParaDiffusionCoef      ChiParaCoef_;
      StateVariableConstantCoef        ChiPerpCoef_;
      StateVariableCoef *              ChiParaCoefPtr_;
      StateVariableCoef *              ChiPerpCoefPtr_;
      Aniso2DDiffusionCoef             ChiCoef_;
      ProductCoefficient               nChiParaCoef_;
      ProductCoefficient               nChiPerpCoef_;
      ProductCoefficient               nkChiParaCoef_;
      ProductCoefficient               nkChiPerpCoef_;
      StateVariableScalarMatrixProductCoef nChiCoef_;
      StateVariableScalarMatrixProductCoef nkChiCoef_;
      StateVariableScalarVectorProductCoef keVCoef_;

      IonEnergyIonizationCoef     SIZCoef_;
      IonEnergyRecombinationCoef  SRCCoef_;
      IonEnergyChargeExchangeCoef SCXCoef_;

      IonAdvectionCoef            ViCoef_;

      ParGridFunction * ChiParaGF_;
      ParGridFunction * ChiPerpGF_;
      ParGridFunction * AdvGF_;
      ParGridFunction * SIZGF_;
      ParGridFunction * SRCGF_;
      ParGridFunction * SCXGF_;
      ParGridFunction * SGF_;
      ParGridFunction * QiGF_;
      ParGridFunction * totEnergyGF_;

   public:
      friend class IonTotalEnergyTerms;

      IonTotalEnergyOp(const DGParams & dg,
                       const PlasmaParams & plasma,
                       ParFiniteElementSpace & vfes,
                       ParFiniteElementSpace & h1_fes,
                       ParGridFunctionArray & yGF,
                       ParGridFunctionArray & kGF,
                       ParGridFunction & elOrdGF,
                       ParGridFunction & hGF,
                       const AdvectionDiffusionBC & bcs,
                       const CoupledBCs & cbcs,
                       const ITECoefs & espcoefs,
                       const CmnCoefs & cmncoefs,
                       VectorCoefficient & B3Coef,
                       double ChiPerp,
                       int term_flag, int vis_flag,
                       int logging,
                       const std::string & log_prefix);

      virtual ~IonTotalEnergyOp();

      virtual void SetTime(double t);

      void Update();

      virtual int GetDefaultVisFlag() { return 255; }

      virtual void RegisterDataFields(DataCollection & dc);

      virtual void PrepareDataFields();
   };

   class ElectronTotalEnergyOp : public TotalEnergyOp
   {
   private:
      enum TermFlag {DIFFUSION_TERM = 0,
                     ADVECTION_TERM, KE_ADVECTION_TERM,
                     IONIZATION_SINK_TERM, RECOMBINATION_SINK_TERM,
                     EQUIPARTITION_SOURCE_TERM, SOURCE_TERM, NUM_TERMS
                    };

      enum VisField {DIFFUSION_PARA_COEF = 1, DIFFUSION_PERP_COEF,
                     ADVECTION_COEF, IONIZATION_SINK_COEF,
                     RECOMBINATION_SINK_COEF, EQUIPARTITION_SOURCE_COEF,
                     SOURCE_COEF, ELECTRON_TOTAL_ENERGY, NUM_FIELDS_PLUS_ONE
                    };

      const ElectronTotalEnergyCoefs & etecoefs_;

      double ChiPerpConst_;

      TotalEnergyCoef                  totEnergyCoef_;
      KineticEnergyCoef                kinEnergyCoef_;
      TotalEnergyAdvectionCoef         advFluxCoef_;
      StaticPressureAdvectionCoef      aneViCoef_;
      ElectronThermalParaDiffusionCoef ChiParaCoef_;
      StateVariableConstantCoef        ChiPerpCoef_;
      StateVariableCoef *              ChiParaCoefPtr_;
      StateVariableCoef *              ChiPerpCoefPtr_;
      Aniso2DDiffusionCoef             ChiCoef_;
      ProductCoefficient               nChiParaCoef_;
      ProductCoefficient               nChiPerpCoef_;
      ProductCoefficient               nkChiParaCoef_;
      ProductCoefficient               nkChiPerpCoef_;
      StateVariableScalarMatrixProductCoef nChiCoef_;
      StateVariableScalarMatrixProductCoef nkChiCoef_;
      StateVariableScalarVectorProductCoef keVCoef_;

      ElectronEnergyIonizationCoef     SIZCoef_;
      ElectronEnergyRecombinationCoef  SRCCoef_;

      ParGridFunction * ChiParaGF_;
      ParGridFunction * ChiPerpGF_;
      ParGridFunction * AdvGF_;
      ParGridFunction * SIZGF_;
      ParGridFunction * SRCGF_;
      ParGridFunction * SGF_;
      ParGridFunction * QiGF_;
      ParGridFunction * totEnergyGF_;

   public:
      friend class ElectronTotalEnergyTerms;

      ElectronTotalEnergyOp(const DGParams & dg,
                            const PlasmaParams & plasma,
                            ParFiniteElementSpace & vfes,
                            ParFiniteElementSpace & h1_fes,
                            ParGridFunctionArray & yGF,
                            ParGridFunctionArray & kGF,
                            ParGridFunction & elOrdGF,
                            ParGridFunction & hGF,
                            const AdvectionDiffusionBC & bcs,
                            const CoupledBCs & cbcs,
                            const ETECoefs & espcoefs,
                            const CmnCoefs & cmncoefs,
                            VectorCoefficient & B3Coef,
                            double ChiPerp,
                            int term_flag, int vis_flag,
                            int logging,
                            const std::string & log_prefix);

      virtual ~ElectronTotalEnergyOp();

      virtual void SetTime(double t);

      void Update();

      virtual int GetDefaultVisFlag() { return 15; }

      virtual void RegisterDataFields(DataCollection & dc);

      virtual void PrepareDataFields();
   };

   class DummyOp : public TransportOp
   {
   public:
      DummyOp(const DGParams & dg,
              const PlasmaParams & plasma,
              ParGridFunctionArray & yGF,
              ParGridFunctionArray & kGF,
              ParGridFunction & elOrdGF,
              ParGridFunction & hGF,
              const AdvectionDiffusionBC & bcs,
              const CoupledBCs & cbcs,
              const CommonCoefs & cmncoefs,
              VectorCoefficient & B3Coef,
              int index,
              const std::string &eqn_name,
              const std::string &field_name,
              int term_flag, int vis_flag,
              int logging, const std::string & log_prefix);

      virtual void SetTimeStep(double dt)
      {
         if (Mpi::Root() && logging_ > 1)
         {
            std::cout << "Setting time step: " << dt << " in DummyOp\n";
         }
         TransportOp::SetTimeStep(dt);
      }

      void Update();

      virtual int GetDefaultVisFlag() { return 0; }
   };

   class VisualizationOp;

   class CombinedOp : public Operator
   {
   private:
      int neq_;
      int logging_;
      int op_flag_;

      ParFiniteElementSpace &fes_;
      // ParGridFunctionArray  &yGF_;
      ParGridFunctionArray  &kGF_;

      Array<TransportOp*> op_;

      const Vector &wgts_;

      Array<int> & offsets_;
      mutable BlockOperator *grad_;

      common::L2_ParFESpace   l2_fes_0_;
      common::H1_ParFESpace   h1_fes_1_;

      VectorCoefficient & B3Coef_;

      mutable ParGridFunction         elOrdDiscGF_;
      mutable GridFunctionCoefficient elOrdDiscCoef_;
      mutable ParGridFunction         elOrdContGF_;

      mutable ParGridFunction         hDiscGF_;
      mutable GridFunctionCoefficient hDiscCoef_;
      mutable ParGridFunction         hContGF_;

      void updateOffsets();

   public:
      CombinedOp(const DGParams & dg,
                 const std::vector<ArtViscParams> & av,
                 const PlasmaParams & plasma, const Vector &eqn_weights,
                 ParFiniteElementSpace & vfes,
                 ParFiniteElementSpace & h1_fes,
                 ParGridFunctionArray & yGF, ParGridFunctionArray & kGF,
                 const TransportBCs & bcs,
                 const TransportCoefs & coefs,
                 Array<int> & offsets,
                 double DiPerp, double XiPerp, double XePerp,
                 const Array<int> & term_flags,
                 const Array<int> & vis_flags,
                 unsigned int op_flag = 31, int logging = 0);

      ~CombinedOp();

      void SetTime(double t);
      void SetTimeStep(double dt);
      void SetLogging(int logging);

      bool IsEquationActive(int eq) const { return  (op_flag_ >> eq) & 1; }

      inline Coefficient * GetDnCoef()
      { return op_[0]->GetDiffusionCoef(); }
      inline MatrixCoefficient * GetDiCoef()
      { return op_[1]->GetDiffusionMatrixCoef(); }
      inline MatrixCoefficient * GetEtaCoef()
      { return op_[2]->GetDiffusionMatrixCoef(); }
      inline MatrixCoefficient * GetnXiCoef()
      { return op_[3]->GetDiffusionMatrixCoef(); }
      inline MatrixCoefficient * GetnXeCoef()
      { return op_[4]->GetDiffusionMatrixCoef(); }

      inline Solver* GetPreconditionerBlock(const TransPrecParams &pparams,
                                            int i) const
      {
         return op_[i]->GetPreconditioner(pparams);
      }

      void Update();

      void Mult(const Vector &k, Vector &y) const;

      void UpdateGradient(const Vector &x) const;

      Operator &GetGradient(const Vector &x) const
      { UpdateGradient(x); return *grad_; }

      void RegisterDataFields(DataCollection & dc);

      void PrepareDataFields();

      void InitializeGLVis();

      void DisplayToGLVis();
   };

   class TransportLeftPrec : public BlockDiagonalPreconditioner
   {
   private:
      int logging_;

      TransPrecParams p_;

      Array<Operator*> diag_prec_;
#ifdef MFEM_USE_SUPERLU
      Array<SuperLURowLocMatrix*> slu_mat_;
#endif
#ifdef MFEM_USE_STRUMPACK
      Array<STRUMPACKRowLocMatrix*> stp_mat_;
#endif
      DGTransportTDO::CombinedOp & comb_op_;

   public:
      TransportLeftPrec(const TransPrecParams &p,
                        const Array<int> &offsets,
                        CombinedOp &combOp,
                        int logging = 0);
      ~TransportLeftPrec();

      virtual void SetOperator(const Operator &op);
   };

   class RightPreconditioner
   {
   protected:
      int height; ///< Dimension of the output / number of rows in the matrix.
      int width;  ///< Dimension of the input / number of columns in the matrix.

   public:
      RightPreconditioner(int size) : height(size), width(size) {}

      /// Get the height (size of output) of the Operator. Synonym with NumRows().
      inline int Height() const { return height; }
      /** @brief Get the number of rows (size of output) of the Operator. Synonym
      with Height(). */
      inline int NumRows() const { return height; }

      /// Get the width (size of input) of the Operator. Synonym with NumCols().
      inline int Width() const { return width; }
      /** @brief Get the number of columns (size of input) of the Operator. Synonym
      with Width(). */
      inline int NumCols() const { return width; }

      virtual void Mult(const Vector&x, Vector &y) = 0;
      virtual void InverseMult(const Vector&y, Vector &x) = 0;
   };

   class RightBlockDiagonalPreconditioner : public RightPreconditioner
   {
   protected:
      Array<int> offsets_;
      const Vector &scale_factors_;

   public:
      RightBlockDiagonalPreconditioner(const Array<int> &offsets,
                                       const Vector &scale_factors)
         : RightPreconditioner(offsets[offsets.Size()-1]),
           offsets_(offsets), scale_factors_(scale_factors)
      {
         MFEM_VERIFY(offsets.Size() - 1 == scale_factors.Size(),
                     "RightBlockDiagonalPreconditioner: "
                     "Incompatible numbers of offsets and scale factors.");
      }

      //! Return the offsets for block starts
      Array<int> & Offsets() { return offsets_; }

      //! Read only access to the offsets for block starts
      const Array<int> & Offsets() const { return offsets_; }

      virtual void Mult(const Vector&x, Vector &y)
      {
         for (int i=0; i<offsets_.Size() - 1; i++)
         {
            for (int j=offsets_[i]; j<offsets_[i+1]; j++)
            {
               y[j] = scale_factors_[i] * x[j];
            }
         }
      }

      virtual void InverseMult(const Vector&y, Vector &x)
      {
         for (int i=0; i<offsets_.Size() - 1; i++)
         {
            for (int j=offsets_[i]; j<offsets_[i+1]; j++)
            {
               x[j] = y[j] / scale_factors_[i];
            }
         }
      }
   };

   class TransportRightPrec : public RightBlockDiagonalPreconditioner
   {
   private:
      int logging_;

      TransPrecParams p_;

      // DGTransportTDO::CombinedOp & comb_op_;

   public:
      TransportRightPrec(const TransPrecParams &p,
                         const Array<int> &offsets,
                         const Vector &scale_factors,
                         CombinedOp &combOp,
                         int logging = 0)
         : RightBlockDiagonalPreconditioner(offsets, scale_factors),
           logging_(logging), p_(p)
      {}

      ~TransportRightPrec() {}

      virtual void SetOperator(const Operator &op);

      virtual void Mult(const Vector&x, Vector &y);
      virtual void InverseMult(const Vector&y, Vector &x);
   };

   /// GMRES method with right preconditioner
   class GMRESRPCSolver : public IterativeSolver
   {
   protected:
      int m; // see SetKDim()

      RightPreconditioner *r_prec;

      static inline void GeneratePlaneRotation(double &dx, double &dy,
                                               double &cs, double &sn);

      static inline void ApplyPlaneRotation(double &dx, double &dy,
                                            double &cs, double &sn);

      static inline void Update(Vector &x, int k, DenseMatrix &h, Vector &s,
                                Array<Vector*> &v);

   public:
      GMRESRPCSolver() { m = 50; r_prec = NULL; }

#ifdef MFEM_USE_MPI
      GMRESRPCSolver(MPI_Comm comm_)
         : IterativeSolver(comm_) { m = 50; r_prec = NULL; }
#endif

      /// Set the number of iteration to perform between restarts, default is 50.
      void SetKDim(int dim) { m = dim; }

      virtual void SetRightPreconditioner(RightPreconditioner &pr)
      { r_prec = &pr; }

      virtual void Mult(const Vector &b, Vector &x) const;
   };

   CombinedOp op_;

   TransportLeftPrec  newton_op_l_prec_;
   TransportRightPrec newton_op_r_prec_;
   GMRESRPCSolver     newton_op_solver_;
   NewtonSolver       newton_solver_;

   mutable Vector x_;
   mutable Vector y_;
   Vector u_;
   Vector dudt_;

public:
   friend class NeutralDensityTerms;
   friend class IonDensityTerms;
   friend class IonMomentumTerms;
   friend class IonTotalEnergyTerms;
   friend class ElectronTotalEnergyTerms;

   DGTransportTDO(const DGParams & dg,
                  const std::vector<ArtViscParams> & av,
                  const PlasmaParams & plasma,
                  const SolverParams & tol,
                  const Vector &eqn_weights,
                  const Vector &fld_weights,
                  ParFiniteElementSpace &fes,
                  ParFiniteElementSpace &vfes,
                  ParFiniteElementSpace &ffes,
                  ParFiniteElementSpace &h1_fes,
                  Array<int> &offsets,
                  ParGridFunctionArray &yGF,
                  ParGridFunctionArray &kGF,
                  const TransportBCs & bcs,
                  const TransportCoefs & coefs,
                  double Di_perp, double Xi_perp, double Xe_perp,
                  const Array<int> & term_flags,
                  const Array<int> & vis_flags,
                  bool imex = true,
                  unsigned int op_flag = 31,
                  int logging = 0);

   ~DGTransportTDO();

   void SetTime(const double _t);
   void SetLogging(int logging);

   bool IsEquationActive(int eq) const { return  (op_flag_ >> eq) & 1; }

   double CheckGradient();

   bool CheckForSteadyState();

   void RegisterDataFields(DataCollection & dc);

   void PrepareDataFields();

   void InitializeGLVis();

   void DisplayToGLVis();

   inline Coefficient * GetDnCoefficient() { return op_.GetDnCoef(); }
   inline MatrixCoefficient * GetDiCoefficient() { return op_.GetDiCoef(); }
   inline MatrixCoefficient * GetEtaCoefficient() { return op_.GetEtaCoef(); }
   inline MatrixCoefficient * GetnXiCoefficient() { return op_.GetnXiCoef(); }
   inline MatrixCoefficient * GetnXeCoefficient() { return op_.GetnXeCoef(); }

   virtual void ImplicitSolve(const double dt, const Vector &y, Vector &k);

   void Update();
};

class DGTransportTDO::VisualizationOp : public AdvTransportOp
{
private:
   enum VisField {B_POLOIDAL = 0,
                  B_TOROIDAL,
                  COULOMB_LOG,
                  TAU_I,
                  TAU_E,
                  IONIZATION_RATE,
                  RECOMBINATION_RATE,
                  CHARGE_EXCHANGE_RATE,
                  ION_SOUND_SPEED
                 };

   // VectorCoefficient & B3Coef_;
   VectorXYCoefficient BxyCoef_;
   VectorZCoefficient  BzCoef_;

   CoulombLogEICoef          lnLambdaCoef_;
   IonCollisionTimeCoef      TauICoef_;
   ElectronCollisionTimeCoef TauECoef_;
   ApproxIonizationRate      SigmaIZCoef_;
   ApproxRecombinationRate   SigmaRCCoef_;
   ApproxChargeExchangeRate  SigmaCXCoef_;
   // SoundSpeedCoef            CsCoef_;

   ParGridFunction * BxyGF_;
   ParGridFunction * BzGF_;
   ParGridFunction * lnLambdaGF_;
   ParGridFunction * TauIGF_;
   ParGridFunction * TauEGF_;
   ParGridFunction * SigmaIZGF_;
   ParGridFunction * SigmaRCGF_;
   ParGridFunction * SigmaCXGF_;
   ParGridFunction * CsGF_;

public:
   VisualizationOp(const DGParams & dg,
                   const PlasmaParams & plasma,
                   ParFiniteElementSpace &vfes,
                   ParGridFunctionArray & yGF,
                   ParGridFunctionArray & kGF,
                   ParGridFunction & elOrdGF,
                   ParGridFunction & hGF,
                   const AdvectionDiffusionBC & bcs,
                   const CoupledBCs & cbcs,
                   const CommonCoefs & cmncoefs,
                   VectorCoefficient & B3Coef,
                   int vis_flag,
                   int logging,
                   const std::string & log_prefix);

   ~VisualizationOp();

   virtual void SetTimeStep(double dt)
   {
      if (Mpi::Root() && logging_ > 1)
      {
         std::cout << "Setting time step: " << dt << " in VisualizationOp\n";
      }
      TransportOp::SetTimeStep(dt);
   }

   void Update();

   inline bool CheckVisFlag(int flag) { return (vis_flag_>> flag) & 1; }

   virtual int GetDefaultVisFlag() { return 3; }

   virtual void RegisterDataFields(DataCollection & dc);

   virtual void PrepareDataFields();
};

class EquationTerms
{
protected:

   bool enabled_;
   unsigned int mask_;
   std::string eqnName_;
   std::vector<std::string> termNames_;

   static bool checkFlag(unsigned int flag, unsigned int mask)
   { return (mask >> flag) & 1; }

   static void   setFlag(unsigned int flag, unsigned int &mask)
   { mask |= (1 << flag); }

   static void clearFlag(unsigned int flag, unsigned int &mask)
   { mask &= ~(1 << flag); }

   virtual void ReadTerms(std::istream &input);

public:
   EquationTerms(const std::string &eqnName, int nTerms)
      : enabled_(false), mask_(0), eqnName_(eqnName), termNames_(nTerms) {}

   virtual ~EquationTerms() {}

   void LoadTerms(std::istream &input)
   { ReadTerms(input); }

   inline bool IsEnabled() const { return enabled_; }
   inline unsigned int GetTermMask() const { return mask_; }
};

class NeutralDensityTerms : public EquationTerms
{
private:
   typedef DGTransportTDO::NeutralDensityOp OP;

public:
   NeutralDensityTerms();
};

class IonDensityTerms : public EquationTerms
{
private:
   typedef DGTransportTDO::IonDensityOp OP;

public:
   IonDensityTerms();
};

class IonMomentumTerms : public EquationTerms
{
private:
   typedef DGTransportTDO::IonMomentumOp OP;

public:
   IonMomentumTerms();
};

class IonTotalEnergyTerms : public EquationTerms
{
private:
   typedef DGTransportTDO::IonTotalEnergyOp OP;

public:
   IonTotalEnergyTerms();
};

class ElectronTotalEnergyTerms : public EquationTerms
{
private:
   typedef DGTransportTDO::ElectronTotalEnergyOp OP;

public:
   ElectronTotalEnergyTerms();
};

class TransportEquationTerms
{
private:
   int neqn_;
   int logging_;
   Array<EquationTerms *> eqnTerms_;

   int eqnMask_;
   Array<int> eqnTermMasks_;

   static void   setFlag(unsigned int flag, int &mask)
   { mask |= (1 << flag); }

   static void clearFlag(unsigned int flag, int &mask)
   { mask &= ~(1 << flag); }

   void Read(std::istream &input);

public:
   TransportEquationTerms(int neqn, int logging = 0);

   void LoadEquationTerms(std::istream &input)
   { Read(input); }

   int & GetEquationBitMask() { return eqnMask_; }

   Array<int> & GetTermFlags() { return eqnTermMasks_; }

   int& operator()(int i) { return eqnTermMasks_[i]; }
   const int operator()(int i) const { return eqnTermMasks_[i]; }

   int& operator[](int i) { return eqnTermMasks_[i]; }
   const int operator[](int i) const { return eqnTermMasks_[i]; }

};

class MultiSpeciesDiffusion;
class MultiSpeciesAdvection;
/*
class TransportSolver : public ODESolver
{
private:
   ODESolver * impSolver_;
   ODESolver * expSolver_;

   ParFiniteElementSpace & sfes_; // Scalar fields
   ParFiniteElementSpace & vfes_; // Vector fields
   ParFiniteElementSpace & ffes_; // Full system

   BlockVector & nBV_;

   ParGridFunction & B_;

   Array<int> & charges_;
   Vector & masses_;

   MultiSpeciesDiffusion * msDiff_;

   void initDiffusion();

public:
   TransportSolver(ODESolver * implicitSolver, ODESolver * explicitSolver,
                   ParFiniteElementSpace & sfes,
                   ParFiniteElementSpace & vfes,
                   ParFiniteElementSpace & ffes,
                   BlockVector & nBV,
                   ParGridFunction & B,
                   Array<int> & charges,
                   Vector & masses);
   ~TransportSolver();

   void Update();

   void Step(Vector &x, double &t, double &dt);
};
*/
/*
class ChiParaCoefficient : public Coefficient
{
private:
 BlockVector & nBV_;
 ParGridFunction nGF_;
 GridFunctionCoefficient nCoef_;
 GridFunctionCoefficient TCoef_;

 int ion_;
 Array<int> & z_;
 Vector     * m_;
 Vector       n_;

public:
 ChiParaCoefficient(BlockVector & nBV, Array<int> & charges);
 ChiParaCoefficient(BlockVector & nBV, int ion_species,
                    Array<int> & charges, Vector & masses);
 void SetT(ParGridFunction & T);

 double Eval(ElementTransformation &T, const IntegrationPoint &ip);
};

class ChiPerpCoefficient : public Coefficient
{
private:
 int ion_;

public:
 ChiPerpCoefficient(BlockVector & nBV, Array<int> & charges);
 ChiPerpCoefficient(BlockVector & nBV, int ion_species,
                    Array<int> & charges, Vector & masses);

 double Eval(ElementTransformation &T, const IntegrationPoint &ip);
};

class ChiCrossCoefficient : public Coefficient
{
private:
 int ion_;

public:
 ChiCrossCoefficient(BlockVector & nBV, Array<int> & charges);
 ChiCrossCoefficient(BlockVector & nBV, int ion_species,
                     Array<int> & charges, Vector & masses);

 double Eval(ElementTransformation &T, const IntegrationPoint &ip);
};

class ChiCoefficient : public MatrixCoefficient
{
private:
 ChiParaCoefficient  chiParaCoef_;
 ChiPerpCoefficient  chiPerpCoef_;
 ChiCrossCoefficient chiCrossCoef_;
 VectorGridFunctionCoefficient BCoef_;

 Vector bHat_;

public:
 ChiCoefficient(int dim, BlockVector & nBV, Array<int> & charges);
 ChiCoefficient(int dim, BlockVector & nBV, int ion_species,
                Array<int> & charges, Vector & masses);

 void SetT(ParGridFunction & T);
 void SetB(ParGridFunction & B);

 void Eval(DenseMatrix &K, ElementTransformation &T,
           const IntegrationPoint &ip);
};
*/
/*
class EtaParaCoefficient : public Coefficient
{
private:
 BlockVector & nBV_;
 ParGridFunction nGF_;
 GridFunctionCoefficient nCoef_;
 GridFunctionCoefficient TCoef_;

 int ion_;
 Array<int> & z_;
 Vector     * m_;
 Vector       n_;

public:
 EtaParaCoefficient(BlockVector & nBV, Array<int> & charges);
 EtaParaCoefficient(BlockVector & nBV, int ion_species,
                    Array<int> & charges, Vector & masses);

 void SetT(ParGridFunction & T);

 double Eval(ElementTransformation &T, const IntegrationPoint &ip);
};
*/
/*
class MultiSpeciesDiffusion : public TimeDependentOperator
{
private:
   ParFiniteElementSpace &sfes_;
   ParFiniteElementSpace &vfes_;

   BlockVector & nBV_;

   Array<int> & charges_;
   Vector & masses_;

   void initCoefficients();
   void initBilinearForms();

public:
   MultiSpeciesDiffusion(ParFiniteElementSpace & sfes,
                         ParFiniteElementSpace & vfes,
                         BlockVector & nBV,
                         Array<int> & charges,
                         Vector & masses);

   ~MultiSpeciesDiffusion();

   void Assemble();

   void Update();

   void ImplicitSolve(const double dt, const Vector &x, Vector &y);
};
*/

// Time-dependent operator for the right-hand side of the ODE representing the
// DG weak form for the diffusion term. (modified from ex14p)
class DiffusionTDO : public TimeDependentOperator
{
private:
   const int dim_;
   double dt_;
   double dg_sigma_;
   double dg_kappa_;

   ParFiniteElementSpace &fes_;
   // ParFiniteElementSpace &dfes_;
   ParFiniteElementSpace &vfes_;

   ParBilinearForm m_;
   ParBilinearForm d_;

   ParLinearForm rhs_;
   ParGridFunction x_;

   HypreParMatrix * M_;
   HypreParMatrix * D_;

   Vector RHS_;
   Vector X_;

   HypreSolver * solver_;
   HypreSolver * amg_;

   MatrixCoefficient &nuCoef_;
   ScalarMatrixProductCoefficient dtNuCoef_;

   void initSolver(double dt);

public:
   DiffusionTDO(ParFiniteElementSpace &fes,
                ParFiniteElementSpace &dfes,
                ParFiniteElementSpace &_vfes,
                MatrixCoefficient & nuCoef,
                double dg_sigma,
                double dg_kappa);

   // virtual void Mult(const Vector &x, Vector &y) const;

   virtual void ImplicitSolve(const double dt, const Vector &x, Vector &y);

   virtual ~DiffusionTDO() { }
};

// Time-dependent operator for the right-hand side of the ODE representing the
// DG weak form for the advection term.
class AdvectionTDO : public TimeDependentOperator
{
private:
   const int dim_;
   const int num_equation_;
   const double specific_heat_ratio_;

   mutable double max_char_speed_;

   ParFiniteElementSpace &vfes_;
   Operator &A_;
   SparseMatrix &Aflux_;
   DenseTensor Me_inv_;

   mutable Vector state_;
   mutable DenseMatrix f_;
   mutable DenseTensor flux_;
   mutable Vector z_;

   void GetFlux(const DenseMatrix &state, DenseTensor &flux) const;

public:
   AdvectionTDO(ParFiniteElementSpace &_vfes,
                Operator &A, SparseMatrix &Aflux, int num_equation,
                double specific_heat_ratio);

   virtual void Mult(const Vector &x, Vector &y) const;

   virtual ~AdvectionTDO() { }
};

// Implements a simple Rusanov flux
class RiemannSolver
{
private:
   int num_equation_;
   double specific_heat_ratio_;
   Vector flux1_;
   Vector flux2_;

public:
   RiemannSolver(int num_equation, double specific_heat_ratio);
   double Eval(const Vector &state1, const Vector &state2,
               const Vector &nor, Vector &flux);
};


// Constant (in time) mixed bilinear form multiplying the flux grid function.
// The form is (vec(v), grad(w)) where the trial space = vector L2 space (mesh
// dim) and test space = scalar L2 space.
class DomainIntegrator : public BilinearFormIntegrator
{
private:
   Vector shape_;
   DenseMatrix flux_;
   DenseMatrix dshapedr_;
   DenseMatrix dshapedx_;

public:
   DomainIntegrator(const int dim, const int num_equation);

   virtual void AssembleElementMatrix2(const FiniteElement &trial_fe,
                                       const FiniteElement &test_fe,
                                       ElementTransformation &Tr,
                                       DenseMatrix &elmat);
};

// Interior face term: <F.n(u),[w]>
class FaceIntegrator : public NonlinearFormIntegrator
{
private:
   int num_equation_;
   double max_char_speed_;
   RiemannSolver rsolver_;
   Vector shape1_;
   Vector shape2_;
   Vector funval1_;
   Vector funval2_;
   Vector nor_;
   Vector fluxN_;
   IntegrationPoint eip1_;
   IntegrationPoint eip2_;

public:
   FaceIntegrator(RiemannSolver &rsolver_, const int dim,
                  const int num_equation);

   virtual void AssembleFaceVector(const FiniteElement &el1,
                                   const FiniteElement &el2,
                                   FaceElementTransformations &Tr,
                                   const Vector &elfun, Vector &elvect);
};

/** Integrator for the DG form:

    - < {(Q grad(u)).n}, [v] > + sigma < [u], {(Q grad(v)).n} >
    + kappa < {h^{-1} cot(theta_T) Q_max^2/Q_min} [u], [v] >,

    where Q is a matrix diffusion coefficient and u, v are the trial
    and test spaces, respectively. Q_max and Q_min are the global
    approximations of the maximum and minimum eigenvalues of Q. The
    function cot(theta_T) is a measure of the distortion of the mesh
    with theta_T approximating the minimum interior angle of each
    element. The parameters sigma and kappa determine the DG method to
    be used (when this integrator is added to the "broken"
    DiffusionIntegrator):
    * sigma = -1, kappa >= kappa0: symm. interior penalty (IP or SIPG) method,
    * sigma = +1, kappa > 0: non-symmetric interior penalty (NIPG) method,
    * sigma = +1, kappa = 0: the method of Baumann and Oden. */
/*
class DGAnisoDiffusionIntegrator : public BilinearFormIntegrator
{
protected:
   MatrixCoefficient *MQ;
   Coefficient *CotTheta;
   double q0, q1;
   double sigma, kappa;

   // these are not thread-safe!
   Vector shape1, shape2, dshape1dn, dshape2dn, nor, nh, ni;
   DenseMatrix jmat, dshape1, dshape2, mq, adjJ;

public:
   DGAnisoDiffusionIntegrator(MatrixCoefficient &q,
                              Coefficient &cotTheta,
                              const double qMin, const double qMax,
                              const double s, const double k)
      : MQ(&q), CotTheta(&cotTheta), q0(qMin), q1(qMax), sigma(s), kappa(k) { }
   using BilinearFormIntegrator::AssembleFaceMatrix;
   virtual void AssembleFaceMatrix(const FiniteElement &el1,
                                   const FiniteElement &el2,
                                   FaceElementTransformations &Trans,
                                   DenseMatrix &elmat);
};
*/

class DGAdvDiffBaseIntegrator
{
protected:
   Coefficient *Q;
   MatrixCoefficient *MQ;
   VectorCoefficient *beta;
   Coefficient *QPara;
   Coefficient *QPerp;
   double lambda, sigma, kappa1, kappa2;

   DGAdvDiffBaseIntegrator(Coefficient & q, VectorCoefficient & b,
                           double l, double s, double k1, double k2)
      :
      Q(&q),
      MQ(NULL),
      beta(&b),
      QPara(NULL),
      QPerp(NULL),
      lambda(l),
      sigma(s),
      kappa1(k1),
      kappa2(k2)
   { }

   DGAdvDiffBaseIntegrator(MatrixCoefficient & q, VectorCoefficient & b,
                           Coefficient *qPara, Coefficient *qPerp,
                           double l, double s, double k1, double k2)
      :
      Q(NULL),
      MQ(&q),
      beta(&b),
      QPara(qPara),
      QPerp(qPerp),
      lambda(l),
      sigma(s),
      kappa1(k1),
      kappa2(k2)
   { }

};

/** Integrator for the DG form:

    < {- Q grad(u) + beta u}_alpha, [v] >
    - sigma < [u], {- Q grad(v) + beta v}_alpha >
    + sigma < [u], {beta v} >
    + kappa < [u], [v] >

    Where:
       {Psi}_alpha = alpha_1 Psi_1 + alpha_2 Psi_2 and alpha_2 = 1 - alpha_1
       {Psi} = (Psi_1 + Psi_2) / 2
       [phi] = n_1 phi_1 + n_2 phi_2
    The diffusion coefficient is the matrix Q, the advection coefficient is
    the vector beta.  The parameter sigma determines the DG method to be used
    (when this integrator is added to the "broken" DiffusionIntegrator and
    the ConservativeConvectionIntegrator):
    * sigma = -1: symm. interior penalty (IP or SIPG) method,
    * sigma = +1: non-symmetric interior penalty (NIPG) method,

    The alpha parameters are determined using a continuous scalar field tau
    according to:
       alpha = (0.5, 0.5) + 0.5 tau (sign(beta.n_1), sign(beta.n_2))
    When tau = 0 this leads to an equal weighting across interelement
    boundaries. When tau = 1 this leads to classical upwinding. Values between
    these extremes can be used to control the degree of upwinding between each
    pair of elements.

    The parameter kappa is a penalty parameter which encourages continuity of
    the solution. See the 2007 paper "Estimation of penalty parameters for
    symmetric interior penalty Galerkin methods" by Y. Epshteyn and B. Riviere
    for advice on selecting kappa. Crudely the rule is:
       kappa > p (p+1) f(Q) g(T_h) in 2D
       kappa > p (p+2) f(Q) g(T_h) in 3D
    Where g(T_h) is a function of the mesh spacing and distortion, and f(Q)
    is (Q_max)^2/Q_min with Q_max and Q_min being the maximum and minimum
    eigenvalues of the diffusion coefficient. It's likely that the advection
    velocity should also contribute to kappa but exactly how is unclear.

    Finally it should be noted that this formulation is borowed from the 2009
    paper "Discontinuous Galerkin methods for advection-diffusion-reaction
    problems" by B. Ayuso and D. Marini where it appears as equation 3.8 (with
    slight modifications). Note in particular that our sigma parameter is the
    negative of the theta parameter from the paper.
*/
class DGAdvDiffIntegrator : public BilinearFormIntegrator,
   DGAdvDiffBaseIntegrator
{
private:
#ifndef MFEM_THREAD_SAFE
   Vector shape1, shape2, nQdshape1, nQdshape2;
   DenseMatrix dshape1, dshape2;
#endif

   double ComputeUpwindingParam(double epsilon, double betaMag);

public:
   DGAdvDiffIntegrator(Coefficient & q, VectorCoefficient & b,
                       double l, double s, double k1, double k2)
      : DGAdvDiffBaseIntegrator(q, b, l, s, k1, k2) {}

   DGAdvDiffIntegrator(MatrixCoefficient & q, VectorCoefficient & b,
                       Coefficient *qPara, Coefficient *qPerp,
                       double l, double s, double k1, double k2)
      : DGAdvDiffBaseIntegrator(q, b, qPara, qPerp, l, s, k1, k2) {}

   using BilinearFormIntegrator::AssembleFaceMatrix;
   virtual void AssembleFaceMatrix(const FiniteElement &el1,
                                   const FiniteElement &el2,
                                   FaceElementTransformations &Trans,
                                   DenseMatrix &elmat);
};

class DGAdvDiffBdrIntegrator : public BilinearFormIntegrator,
   DGAdvDiffBaseIntegrator
{
private:
#ifndef MFEM_THREAD_SAFE
   Vector shape1, nQdshape1;
   DenseMatrix dshape1;
#endif

public:
   DGAdvDiffBdrIntegrator(Coefficient & q, VectorCoefficient & b,
                          double l, double s, double k1, double k2)
      : DGAdvDiffBaseIntegrator(q, b, l, s, k1, k2) {}

   DGAdvDiffBdrIntegrator(MatrixCoefficient & q, VectorCoefficient & b,
                          Coefficient *qPara, Coefficient *qPerp,
                          double l, double s, double k1, double k2)
      : DGAdvDiffBaseIntegrator(q, b, qPara, qPerp, l, s, k1, k2) {}

   using BilinearFormIntegrator::AssembleFaceMatrix;
   virtual void AssembleFaceMatrix(const FiniteElement &el1,
                                   const FiniteElement &el2,
                                   FaceElementTransformations &Trans,
                                   DenseMatrix &elmat);
};

/** Boundary linear integrator for imposing non-zero Dirichlet boundary
    conditions, to be used in conjunction with DGDiffusionIntegrator.
    Specifically, given the Dirichlet data u_D, the linear form assembles the
    following integrals on the boundary:

    sigma < u_D, (Q grad(v)).n > + kappa < {h^{-1} Q} u_D, v >,

    where Q is a scalar or matrix diffusion coefficient and v is the test
    function. The parameters sigma and kappa should be the same as the ones
    used in the DGDiffusionIntegrator. */
class DGAdvDiffDirichletLFIntegrator : public LinearFormIntegrator,
   DGAdvDiffBaseIntegrator
{
private:
#ifndef MFEM_THREAD_SAFE
   Vector shape1, shape2, nQdshape1, nQdshape2;
   DenseMatrix dshape1, dshape2;
#endif

   Coefficient *uD;

public:
   DGAdvDiffDirichletLFIntegrator(Coefficient &u,
                                  Coefficient &q,
                                  VectorCoefficient & b,
                                  double l,
                                  double s,
                                  double k1,
                                  double k2)
      : DGAdvDiffBaseIntegrator(q, b, l, s, k1, k2),
        uD(&u)
   { }

   DGAdvDiffDirichletLFIntegrator(Coefficient &u,
                                  MatrixCoefficient &q,
                                  VectorCoefficient & b,
                                  Coefficient *qPara,
                                  Coefficient *qPerp,
                                  double l,
                                  double s,
                                  double k1,
                                  double k2)
      : DGAdvDiffBaseIntegrator(q, b, qPara, qPerp, l, s, k1, k2),
        uD(&u)
   { }

   virtual void AssembleRHSElementVect(const FiniteElement &el,
                                       ElementTransformation &Tr,
                                       Vector &elvect)
   {
      mfem_error("DGAdvDiffDirichletLFIntegrator::AssembleRHSElementVect");
   }
   virtual void AssembleRHSElementVect(const FiniteElement &el,
                                       FaceElementTransformations &Tr,
                                       Vector &elvect);
};

/** Integrator for the DG form:

    < {- Q grad(u)}_alpha, [v] >
    - sigma < [u], {- Q grad(v)}_alpha >
    + kappa < [u], [v] >

    Where:
       {Psi}_alpha = alpha_1 Psi_1 + alpha_2 Psi_2 and alpha_2 = 1 - alpha_1
       {Psi} = (Psi_1 + Psi_2) / 2
       [phi] = n_1 phi_1 + n_2 phi_2
    The diffusion coefficient is the matrix Q.  The parameter sigma determines
    the DG method to be used
    (when this integrator is added to the "broken" DiffusionIntegrator:
    * sigma = -1: symm. interior penalty (IP or SIPG) method,
    * sigma = +1: non-symmetric interior penalty (NIPG) method,

    The alpha parameters are determined using a continuous scalar field tau
    according to:
       alpha = (0.5, 0.5) + 0.5 tau (sign(beta.n_1), sign(beta.n_2))

    The parameter kappa is a penalty parameter which encourages continuity of
    the solution. See the 2007 paper "Estimation of penalty parameters for
    symmetric interior penalty Galerkin methods" by Y. Epshteyn and B. Riviere
    for advice on selecting kappa. Crudely the rule is:
       kappa > p (p+1) f(Q) g(T_h) in 2D
       kappa > p (p+2) f(Q) g(T_h) in 3D
    Where g(T_h) is a function of the mesh spacing and distortion, and f(Q)
    is (Q_max)^2/Q_min with Q_max and Q_min being the maximum and minimum
    eigenvalues of the diffusion coefficient. It's likely that the advection
    velocity should also contribute to kappa but exactly how is unclear.

    Finally it should be noted that this formulation is borowed from the 2009
    paper "Discontinuous Galerkin methods for advection-diffusion-reaction
    problems" by B. Ayuso and D. Marini where it appears as equation 3.8 (with
    slight modifications). Note in particular that our sigma parameter is the
    negative of the theta parameter from the paper.
*/
class DGAnisoDiffBaseIntegrator
{
protected:
   MatrixCoefficient *MQ;
   Coefficient *QPara;
   Coefficient *QPerp;
   double sigma, kappa;

   DGAnisoDiffBaseIntegrator(MatrixCoefficient & q,
                             Coefficient *qPara, Coefficient *qPerp,
                             double s, double k)
      :
      MQ(&q),
      QPara(qPara),
      QPerp(qPerp),
      sigma(s),
      kappa(k)
   { }

};

class DGAnisoDiffIntegrator : public BilinearFormIntegrator,
   DGAnisoDiffBaseIntegrator
{
private:
#ifndef MFEM_THREAD_SAFE
   Vector shape1, shape2, nQdshape1, nQdshape2;
   DenseMatrix dshape1, dshape2;
#endif

public:
   DGAnisoDiffIntegrator(MatrixCoefficient & q,
                         Coefficient *qPara, Coefficient *qPerp,
                         double s, double k)
      : DGAnisoDiffBaseIntegrator(q, qPara, qPerp, s, k)
   {}

   using BilinearFormIntegrator::AssembleFaceMatrix;
   virtual void AssembleFaceMatrix(const FiniteElement &el1,
                                   const FiniteElement &el2,
                                   FaceElementTransformations &Trans,
                                   DenseMatrix &elmat);
};

class DGAnisoDiffBdrIntegrator : public BilinearFormIntegrator,
   DGAnisoDiffBaseIntegrator
{
private:
#ifndef MFEM_THREAD_SAFE
   Vector shape1, nQdshape1;
   DenseMatrix dshape1;
#endif

public:
   DGAnisoDiffBdrIntegrator(MatrixCoefficient & q,
                            Coefficient *qPara, Coefficient *qPerp,
                            double s, double k)
      : DGAnisoDiffBaseIntegrator(q, qPara, qPerp, s, k) {}

   using BilinearFormIntegrator::AssembleFaceMatrix;
   virtual void AssembleFaceMatrix(const FiniteElement &el1,
                                   const FiniteElement &el2,
                                   FaceElementTransformations &Trans,
                                   DenseMatrix &elmat);
};

/** Boundary linear integrator for imposing non-zero Dirichlet boundary
    conditions, to be used in conjunction with DGDiffusionIntegrator.
    Specifically, given the Dirichlet data u_D, the linear form assembles the
    following integrals on the boundary:

    sigma < u_D, (Q grad(v)).n > + kappa < {h^{-1} Q} u_D, v >,

    where Q is a scalar or matrix diffusion coefficient and v is the test
    function. The parameters sigma and kappa should be the same as the ones
    used in the DGDiffusionIntegrator. */
class DGAnisoDiffDirichletLFIntegrator : public LinearFormIntegrator,
   DGAnisoDiffBaseIntegrator
{
private:
   Coefficient *uD;

   // these are not thread-safe!
   Vector shape, dshape_dn, nor, nh, ni, vb;
   DenseMatrix dshape, mq, adjJ;

public:
   DGAnisoDiffDirichletLFIntegrator(Coefficient &u,
                                    MatrixCoefficient &q,
                                    Coefficient *qPara,
                                    Coefficient *qPerp,
                                    double s,
                                    double k)
      : DGAnisoDiffBaseIntegrator(q, qPara, qPerp, s, k),
        uD(&u)
   { }

   virtual void AssembleRHSElementVect(const FiniteElement &el,
                                       ElementTransformation &Tr,
                                       Vector &elvect)
   {
      mfem_error("DGAnisoDiffDirichletLFIntegrator::AssembleRHSElementVect");
   }
   virtual void AssembleRHSElementVect(const FiniteElement &el,
                                       FaceElementTransformations &Tr,
                                       Vector &elvect);
};

} // namespace transport

} // namespace plasma

} // namespace mfem

#endif // MFEM_USE_MPI

#endif // MFEM_TRANSPORT_SOLVER
