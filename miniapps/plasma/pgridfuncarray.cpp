// Copyright (c) 2010-2021, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#include "pgridfuncarray.hpp"

#ifdef MFEM_USE_MPI

using namespace std;
namespace mfem
{

namespace plasma
{

ParGridFunctionArray::ParGridFunctionArray(int size, ParFiniteElementSpace *pf)
   : Array<ParGridFunction*>(size), owns_data(true)
{
   for (int i=0; i<size; i++)
   {
      data[i] = new ParGridFunction(pf);
   }
}

ParGridFunctionArray::~ParGridFunctionArray()
{
   if (owns_data)
   {
      for (int i=0; i<size; i++)
      {
         delete data[i];
      }
   }
}

void ParGridFunctionArray::ProjectCoefficient(Array<Coefficient*> &coeff)
{
   for (int i=0; i<size; i++)
   {
      if (coeff[i] != NULL)
      {
         data[i]->ProjectCoefficient(*coeff[i]);
      }
      else
      {
         *data[i] = 0.0;
      }
   }
}

void ParGridFunctionArray::Update()
{
   for (int i=0; i<size; i++)
   {
      data[i]->Update();
   }
}

void ParGridFunctionArray::ExchangeFaceNbrData()
{
   for (int i=0; i<size; i++)
   {
      data[i]->ExchangeFaceNbrData();
   }
}

} // namespace plasma

} // namespace mfem

#endif // MFEM_USE_MPI
