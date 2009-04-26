// Copyright (C) 2001-2008 Vivien Mallet
//
// This file is part of Seldon library.
// Seldon library provides matrices and vectors structures for
// linear algebra.
//
// Seldon is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// Seldon is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License (file "license") for more details.
//
// For more information, please see the Seldon home page:
//     http://spacetown.free.fr/lib/seldon/

#ifndef SELDON_FILE_ITERATIVE_CXX

// headers of class Iteration and Preconditioner_Base
#include "Iterative.hxx"

// and all iterative solvers
#include "Cg.cxx"
#include "Cgne.cxx"
#include "Lsqr.cxx"
#include "Cgs.cxx"
#include "BiCg.cxx"
#include "BiCgStab.cxx"
#include "BiCgStabl.cxx"
#include "BiCgcr.cxx"
#include "Gcr.cxx"
#include "CoCg.cxx"
#include "Gmres.cxx"
#include "MinRes.cxx"
#include "Qmr.cxx"
#include "QmrSym.cxx"
#include "QCgs.cxx"
#include "TfQmr.cxx"
#include "Symmlq.cxx"


namespace Seldon
{
  
  //! Default constructor
  Preconditioner_Base::Preconditioner_Base()
  {
  }
  
  
  //! Solves M z = r
  /*!
    Identity preconditioner M = I
  */
  template<class T2, class Storage2, class Allocator2, class Matrix>
  void Preconditioner_Base::Solve(const Matrix& A,
				  const Vector<T2, Storage2, Allocator2> & r,
				  Vector<T2, Storage2, Allocator2> & z)
  {
    Copy(r,z);
  }
  
  
  //! Solves M^t z = r
  /*!
    Identity preconditioner M = I
  */
  template<class T2, class Storage2, class Allocator2, class Matrix>
  void Preconditioner_Base::
  TransSolve(const Matrix& A, const Vector<T2, Storage2, Allocator2> & r,
	     Vector<T2, Storage2, Allocator2> & z)
  {
    Solve(A, r, z);
  }
  
  
  //! Default constructor
  template<class Titer>
  Iteration<Titer>::Iteration()
  {
    tolerance = 1e-6;
    max_iter = 100;
    nb_iter = 0;
    error_code = 0;
    fail_convergence = false;
    print_level = 1;
    init_guess_null = true;
    type_solver = 0; parameter_restart = 10;
    type_preconditioning = 0;
  }
  
  
  //! Constructor with maximum number of iterations and stopping criterion
  template<class Titer>
  Iteration<Titer>::Iteration(int max_iteration, const Titer& tol)
  {
    max_iter = max_iteration;
    tolerance = tol;
    nb_iter = 0;
    error_code = 0;
    fail_convergence = false;
    print_level = 1;
    init_guess_null = true;
    type_solver = 0; parameter_restart = 10;
    type_preconditioning = 0;
  }
  

  //! Copy constructor
  template<class Titer>
  Iteration<Titer>::Iteration(const Iteration<Titer>& outer)
  {
    tolerance = outer.tolerance; facteur_reste = outer.facteur_reste;
    max_iter = outer.max_iter;
    nb_iter = 0;
    print_level = outer.print_level; error_code = outer.error_code;
    init_guess_null = true;
    type_solver = outer.type_solver;
    parameter_restart = outer.parameter_restart;
    type_preconditioning = outer.type_preconditioning;
  }
  
  
  //! Returns the type of solver
  template<class Titer>
  int Iteration<Titer>::GetTypeSolver() const
  {
    return type_solver;
  }
  
  
  //! Returns the restart parameter
  template<class Titer>
  int Iteration<Titer>::GetRestart() const
  {
    return parameter_restart;
  }
  
  
  //! Returns used coefficient to compute relative residual
  template<class Titer>
  Titer Iteration<Titer>::GetFactor() const
  {
    return facteur_reste;
  }
  
  
  //! Returns stopping criterion
  template<class Titer>
  Titer Iteration<Titer>::GetTolerance() const
  {
    return tolerance;
  }
  
  
  //! Returns the number of iterations
  template<class Titer>
  int Iteration<Titer>::GetNumberIteration() const
  {
    return nb_iter;
  }
  
  
  //! Changes the type of solver and preconditioning
  template<class Titer>
  void Iteration<Titer>::SetSolver(int type_resolution,
				   int param_restart, int type_prec)
  {
    type_solver = type_resolution;
    parameter_restart = param_restart;
    type_preconditioning = type_prec;
  }
  
  
  //! Changes the restart parameter
  template<class Titer>
  void Iteration<Titer>::SetRestart(int m)
  {
    parameter_restart = m;
  }
  
  
  //! Changes the stopping criterion
  template<class Titer>
  void Iteration<Titer>::SetTolerance(Titer stopping_criterion)
  {
    tolerance = stopping_criterion;
  }
  
  
  //! Changes the maximum number of iterations
  template<class Titer>
  void Iteration<Titer>::SetMaxNumberIteration(int max_iteration)
  {
    max_iter=max_iteration;
  }

  
  //! Changes the number of iterations
  template<class Titer>
  void Iteration<Titer>::SetNumberIteration(int nb)
  {
    nb_iter = nb;
  }
  
  
  //! Sets to a normal display (residual each 100 iterations)
  template<class Titer>
  void Iteration<Titer>::ShowMessages()
  {
    print_level = 1;
  }
  
  
  //! Sets to a complete display (residual each iteration)
  template<class Titer>
  void Iteration<Titer>::ShowFullHistory()
  {
    print_level = 6;
  }
  
  
  //! Doesn't display any information
  template<class Titer>
  void Iteration<Titer>::HideMessages()
  {
    print_level = 0;
  }
  
  
  //! Initialization with the right hand side
  template<class Titer> template<class T, class Storage, class Allocator>
  int Iteration<Titer>::Init(const Vector<T, Storage, Allocator>& r)
  {
    Titer norme_rhs = Titer(Norm2(r));
    // test of a null right hand side
    if (norme_rhs == Titer(0))
      return -1;
    
    // coefficient used later to compute relative residual
    facteur_reste = Titer(1)/norme_rhs;
    
    // initialization of iterations
    nb_iter = 0;
    return 0; // successful initialization
  }
  
  
  //! Returns true if it is the first iteration
  template<class Titer>
  inline bool Iteration<Titer>::First() const
  {
    if (nb_iter == 0)
      return true;
    
    return false;
  }
  
  
  //! Returns true if the initial guess is null
  template<class Titer>
  inline bool Iteration<Titer>::IsInitGuess_Null() const
  {
    return init_guess_null;
  }
  
  
  //! Returns true if the iterative solver has reached its end
  template<class Titer> template<class T, class Storage, class Allocator>
  inline bool Iteration<Titer>::
  Finished(const Vector<T, Storage, Allocator>& r) const
  {
    // absolute residual
    Titer reste = Norm2(r);
    // computation of relative residual
    reste = facteur_reste*reste;
    
    // displaying residual if required
    if ((print_level >= 1)&&(nb_iter%100 == 0))
      cout<<"Residu at iteration number "<<
	GetNumberIteration()<<"  "<<reste<<endl;
    else if (print_level >= 6)
      cout<<"Residu at iteration number "<<
	GetNumberIteration()<<"  "<<reste<<endl;
    
    // end of iterative solver when residual is small enough
    // or when the number of iterations is too high
    if ((reste < tolerance)||(nb_iter >= max_iter))
      return true;
    
    return false;
  }
  
  
  //! Returns true if the iterative solver has reached its end
  template<class Titer>
  inline bool Iteration<Titer>::Finished(const Titer& r) const
  {
    // relative residual
    Titer reste = facteur_reste*r;
    
    // displaying residual if required
    if ((print_level >= 1)&&(nb_iter%100 == 0))
      cout<<"Residu at iteration number "<<
	GetNumberIteration()<<"  "<<reste<<endl;
    else if (print_level >= 6)
      cout<<"Residu at iteration number "<<
	GetNumberIteration()<<"  "<<reste<<endl;
    
    // end of iterative solver when residual is small enough
    // or when the number of iterations is too high
    if ((reste < tolerance)||(nb_iter >= max_iter))
      return true;
    
    return false;
  }
  
  
  //! Informs of a failure in the iterative solver
  template<class Titer>
  void Iteration<Titer>::Fail(int i, const string& s)
  {
    fail_convergence = true;
    error_code = i;
    // displays information if required
    if ((print_level >= 1)&&(nb_iter%100==0))
      cout<<"Error during resolution : "<<s<<endl;
  }
  
  
  //! Increment the number of iterations
  template<class Titer>
  inline Iteration<Titer>& Iteration<Titer>::operator ++ (void)
  {
    ++nb_iter;
    return *this;
  }
  
  
  //! Returns the error code (if an error occured)
  template<class Titer>
  int Iteration<Titer>::ErrorCode() const
  {
    if (nb_iter >= max_iter)
      return -2;
    
    if (fail_convergence)
      return error_code;
    
    return 0;
  }
  
} // end namespace

#define SELDON_FILE_ITERATIVE_CXX
#endif

