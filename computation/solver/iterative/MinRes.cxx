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

#ifndef SELDON_FILE_ITERATIVE_MINRES_CXX

namespace Seldon
{
  
  //! Solves a linear system by using Minimum Residual (MinRes)
  /*!
    Solves a real symmetric linear system Ax = b with restarted Preconditioned
    Minimum Residual Algorithm.
       
    return value of 0 indicates convergence within the
    maximum number of iterations (determined by the iter object).
    return value of 1 indicates a failure to converge.
    
    See C. PAIGE AND M. SAUNDERS,
    Solution of sparse indefinite systems of linear equations,
    SIAM J. Numer. Anal., 12 (1975), pp. 617-629.
    
    \param[in] A  Real Symmetric Matrix
    \param[in,out] x  Vector on input it is the initial guess
    on output it is the solution
    \param[in] b  Vector right hand side of the linear system
    \param[in] M Right preconditioner
    \param[in] iter Iteration parameters
  */
  template <class Titer, class Matrix, class Vector, class Preconditioner>
  int MinRes(Matrix& A, Vector& x, const Vector& b,
	     Preconditioner& M, Iteration<Titer> & iter)
  {
    const int N = A.GetM();
    if (N <= 0)
      return 0;
    
    typedef typename Vector::value_type Complexe;
    Vector u_old(N), u(N), r(N), v_old(N), v(N),
      w_old(N), w(N), z(N), w_oold(N);
    
    Complexe dp, beta, ibeta, beta_old, alpha, eta, ceta;
    Complexe cold, coold, c, soold, sold, s, rho0, rho1, rho2, rho3;
    
    int success_init = iter.Init(b);
    if (success_init != 0)
      return iter.ErrorCode();
    
    Seldon::Copy(b,r);
    // r = b - A x
    if (!iter.IsInitGuess_Null())
      MltAdd(Complexe(-1), A, x, Complexe(1), r);
    else
      x.Zero();
    
    u_old.Zero(); v_old.Zero(); w_old.Zero(); w.Zero(); w_oold.Zero();
    // preconditioning
    M.Solve(A, r, z);
    dp = DotProd(r, z);
    dp = sqrt(dp); beta = dp; eta = beta;
    Seldon::Copy(r, v); Seldon::Copy(z, u);
    
    ibeta = 1.0 / beta;
    Mlt(ibeta, v); Mlt(ibeta, u);
    
    c = 1.0; s = 0.0; cold = 1.0; sold = 0.0;
    Titer np = Norm2(b);
    
    iter.SetNumberIteration(0);
    // Loop until the stopping criteria are satisfied
    while (!iter.Finished(np))
      {
	// matrix-vector product r = A*u
	Mlt(A, u, r);
	alpha = DotProd(r, u);
	// preconditioning
	M.Solve(A, r, z);
	
	//  r = r - alpha v
	//  z = z - alpha u
	Seldon::Add(-alpha, v, r);
	Seldon::Add(-alpha, u, z);
	//  r = r - beta v_old
	//  z = z - beta u_old
	Seldon::Add(-beta, v_old, r);
	Seldon::Add(-beta, u_old, z);
	
	beta_old = beta;
	
	dp = DotProd(r, z);
	beta = sqrt(dp);
	
	// QR factorization
	coold = cold; cold = c; soold = sold; sold = s;

	rho0 = cold * alpha - coold * sold * beta_old;
	rho1 = sqrt(rho0*rho0 + beta*beta);
	rho2 = sold * alpha + coold * cold * beta_old;
	rho3 = soold * beta_old;
	 
	// Givens rotation
	if (rho1 == Complexe(0) )
	  {
	    iter.Fail(1, "Minres breakdown #1");
	    break;
	  }
	c = rho0 / rho1;
	s = beta / rho1;
	 
	// update
	Seldon::Copy(w_old, w_oold); Seldon::Copy(w, w_old);
	Seldon::Copy(u, w);
	 
	Seldon::Add(-rho2, w_old, w);
	Seldon::Add(-rho3, w_oold, w);
	Mlt(Complexe(1./rho1), w);
	 
	ceta = c*eta;
	Seldon::Add(ceta, w, x);
	eta = -s*eta;
	 
	Seldon::Copy(v, v_old); Seldon::Copy(u, u_old);
	Seldon::Copy(r, v); Seldon::Copy(z, u);
	if (beta == Complexe(0) )
	  {
	    iter.Fail(2, "MinRes breakdown #2");
	    break;
	  }
	ibeta = 1.0/beta;
	Mlt(ibeta, v); Mlt(ibeta, u);
	 
	// residual norm
	np *= abs(s);
	++ iter;
      }
    return iter.ErrorCode();
  }

} // end namespace
					 
#define SELDON_FILE_ITERATIVE_MINRES_CXX
#endif