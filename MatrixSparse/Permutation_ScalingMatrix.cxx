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

#ifndef SELDON_FILE_PERMUTATION_SCALING_MATRIX_CXX

/*
  Functions defined in this file:
  
  PermuteMatrix(A, I, J)
  
  ScaleMatrix(A, Drow, Dcol)
  
  ScaleLeftMatrix(A, Drow)
*/

namespace Seldon
{
  
  //! Permutation of a general matrix stored by rows
  /*!
    B(row_perm(i),col_perm(j)) = A(i,j)  and A = B
    equivalent Matlab operation : A(row_perm, col_perm) = A
  */
  template<class T, class Allocator>
  void PermuteMatrix(Matrix<T, General, ArrayRowSparse, Allocator>& A,
		     const IVect& row_perm, const IVect& col_perm)
  {
    int m = A.GetM(), n, i, i_, j, i2;
    // Vector<T,Vect_Full,Allocator> tmp;
    IVect ind_tmp, iperm(m), rperm(m);
    for (i = 0; i < m; i++)
      {
	iperm(i) = i;
	rperm(i) = i;
      }
    // A(rperm(i),:) will be the place where is the initial row i
    
    // algorithm avoiding the allocation of another matrix
    for (i = 0; i < m; i++)
      {
	// we get the index of row where the row initially placed on row i is
	i2 = rperm(i);
	// we get the new index of this row
	i_ = row_perm(i);
	
	// we fill ind_tmp of the permuted indices of columns of row i
	n = A.GetRowSize(i2);
	ind_tmp.Reallocate(n);
	for (j = 0; j < n; j++)
	  ind_tmp(j) = col_perm(A.Index(i2,j));
	
	// we swap the two rows i and its destination row_perm(i)
	A.SwapRow(i2, i_);
	A.ReplaceIndexRow(i_, ind_tmp);
	
	// we update the indices iperm and rperm
	// in order to keep in memory the place where the row row_perm(i) is
	int i_tmp = iperm(i_);
	iperm(i_) = iperm(i2);
	iperm(i2) = i_tmp;
	rperm(iperm(i_)) = i_;
	rperm(iperm(i2)) = i2;
	
	// we assemble the row i
	A.AssembleRow(i_);
      }
  }
  
  
  //! Permutation of a symmetric matrix stored by columns
  /*!
    B(row_perm(i),col_perm(j)) = A(i,j)  and A = B
    equivalent Matlab operation : A(row_perm, col_perm) = A
  */
  template<class T,class Allocator>
  void PermuteMatrix(Matrix<T, Symmetric, ArrayColSymSparse, Allocator>& A,
		     const IVect& row_perm, const IVect& col_perm)
  {
    // it is assumed that the permuted matrix is still symmetric !
    // for example, the user can provide row_perm = col_perm
    int m = A.GetM(); int nnz = A.GetDataSize();
    IVect IndRow(nnz), IndCol(nnz); Vector<T, Vect_Full, Allocator> Val(nnz);
    
    // first we convert the matrix in coordinate format
    // and we permute the indices
    // IndRow -> indices of the permuted rows
    // IndCol -> indices of the permuted columns
    int k = 0;
    for (int i = 0; i < m; i++)
      {
	for (int j = 0; j < A.GetColumnSize(i); j++)
	  {
	    IndCol(k) = col_perm(i);
	    Val(k) = A.Value(i,j);
	    IndRow(k) = row_perm(A.Index(i,j));
	    if (IndCol(k) <= IndRow(k))
	      {
		// we store only the superior part of the symmetric matrix
		int ind_tmp = IndRow(k);
		IndRow(k) = IndCol(k);
		IndCol(k) = ind_tmp;
	      }
	    k++;
	  }
      }
    
    // we sort with respect to column numbers
    Sort(nnz, IndCol, IndRow, Val);
    
    // A is filled
    k = 0;
    for (int i = 0; i < m; i++)
      {
	int first_index = k;
	// we get the size of the column i
	while ((k < nnz)&&(IndCol(k) <= i))
	  k++;
	
	int size_column = k-first_index;
	// if column not empty
	if (size_column > 0)
	  {
	    A.ReallocateColumn(i, size_column);
	    k = first_index;
	    Sort(k, k+size_column-1, IndRow, Val);
	    for (int j = 0; j < size_column; j++)
	      {
		A.Index(i,j) = IndRow(k);
		A.Value(i,j) = Val(k);
		k++;
	      }
	  }
	else
	  A.ClearColumn(i);
      }
  }
  
  
  //! Permutation of a symmetric matrix stored by rows
  /*!
    B(row_perm(i),col_perm(j)) = A(i,j)  and A = B
    equivalent Matlab operation : A(row_perm, col_perm) = A
  */
  template<class T, class Allocator>
  void PermuteMatrix(Matrix<T, Symmetric,
		     ArrayRowSymComplexSparse, Allocator>& A,
		     const IVect& row_perm,const IVect& col_perm)
  {
    // it is assumed that the permuted matrix is still symmetric !
    // by example, the user can provide row_perm = col_perm
    int m = A.GetM();
    int nnz_real = A.GetRealDataSize(), nnz_imag = A.GetImagDataSize();
    IVect IndRow(nnz_real), IndCol(nnz_real);
    Vector<T, Vect_Full, Allocator> Val(nnz_real);
    
    // first we convert the matrix in coordinate format
    // and we permute the indices
    // IndRow -> indices of the permuted rows
    // IndCol -> indices of the permuted columns
    int k = 0;
    for (int i = 0; i < m; i++)
      {
	for (int j = 0; j < A.GetRealRowSize(i); j++)
	  {
	    IndRow(k) = row_perm(i);
	    Val(k) = A.ValueReal(i,j);
	    IndCol(k) = col_perm(A.IndexReal(i,j));
	    if (IndCol(k) <= IndRow(k))
	      {
		// we store only the superior part of the symmetric matrix
		int ind_tmp = IndRow(k);
		IndRow(k) = IndCol(k);
		IndCol(k) = ind_tmp;
	      }
	    k++;
	  }
      }
    
    // we sort by row number
    Sort(nnz_real, IndRow, IndCol, Val);
    
    // A is filled
    k = 0;
    for (int i = 0; i < m; i++)
      {
	int first_index = k;
	// we get the size of the row i
	while ((k < nnz_real)&&(IndRow(k) <= i))
	  k++;
	
	int size_row = k - first_index;
	// if row not empty
	if (size_row > 0)
	  {
	    A.ReallocateRealRow(i, size_row);
	    k = first_index;
	    Sort(k, k+size_row-1, IndCol, Val);
	    for (int j = 0; j < size_row; j++)
	      {
		A.IndexReal(i,j) = IndCol(k);
		A.ValueReal(i,j) = Val(k);
		k++;
	      }
	  }
	else
	  A.ClearRealRow(i);
      }
    
    // same procedure for imaginary part

    IndRow.Reallocate(nnz_imag); IndCol.Reallocate(nnz_imag);
    Val.Reallocate(nnz_imag);
    
    // first we convert the matrix in coordinate format
    // and we permute the indices
    // IndRow -> indices of the permuted rows
    // IndCol -> indices of the permuted columns
    k = 0;
    for (int i = 0; i < m; i++)
      {
	for (int j = 0; j < A.GetImagRowSize(i); j++)
	  {
	    IndRow(k) = row_perm(i);
	    Val(k) = A.ValueImag(i,j);
	    IndCol(k) = col_perm(A.IndexImag(i,j));
	    if (IndCol(k) <= IndRow(k))
	      {
		// we store only the superior part of the symmetric matrix
		int ind_tmp = IndRow(k);
		IndRow(k) = IndCol(k);
		IndCol(k) = ind_tmp;
	      }
	    k++;
	  }
      }
    // we sort by row number
    Sort(nnz_imag, IndRow, IndCol, Val);
    
    // A is filled
    k = 0;
    for (int i = 0; i < m; i++)
      {
	int first_index = k;
	// we get the size of the row i
	while ((k < nnz_imag)&&(IndRow(k) <= i))
	  k++;
	int size_row = k - first_index;
	// if row not empty
	if (size_row > 0)
	  {
	    A.ReallocateImagRow(i,size_row);
	    k = first_index;
	    Sort(k,k+size_row-1,IndCol,Val);
	    for (int j = 0; j < size_row; j++)
	      {
		A.IndexImag(i,j) = IndCol(k);
		A.ValueImag(i,j) = Val(k);
		k++;
	      }
	  }
	else
	  A.ClearImagRow(i);
      }
  }
  
  
  //! Permutation of a symmetric matrix stored by rows
  /*!
    B(row_perm(i),col_perm(j)) = A(i,j)  and A = B
    equivalent Matlab operation : A(row_perm, col_perm) = A
  */
  template<class T,class Allocator>
  void PermuteMatrix(Matrix<T, Symmetric, ArrayRowSymSparse, Allocator>& A,
		     const IVect& row_perm, const IVect& col_perm)
  {
    // it is assumed that the permuted matrix is still symmetric !
    // by example, the user can provide row_perm = col_perm
    int m = A.GetM(); int nnz = A.GetDataSize();
    IVect IndRow(nnz),IndCol(nnz); Vector<T,Vect_Full,Allocator> Val(nnz);
    
    // first we convert the matrix in coordinate format
    // and we permute the indices
    // IndRow -> indices of the permuted rows
    // IndCol -> indices of the permuted columns
    int k = 0;
    for (int i = 0; i < m; i++)
      {
	for (int j = 0; j < A.GetRowSize(i); j++)
	  {
	    IndRow(k) = row_perm(i);
	    Val(k) = A.Value(i,j);
	    IndCol(k) = col_perm(A.Index(i,j));
	    if (IndCol(k) <= IndRow(k))
	      {
		// we store only the superior part of the symmetric matrix
		int ind_tmp = IndRow(k);
		IndRow(k) = IndCol(k);
		IndCol(k) = ind_tmp;
	      }
	    k++;
	  }
      }
    // we sort with respect to row numbers
    Sort(nnz, IndRow, IndCol, Val);
    
    // A is filled
    k = 0;
    for (int i = 0; i < m; i++)
      {
	int first_index = k;
	// we get the size of the row i
	while ((k < nnz)&&(IndRow(k) <= i))
	  k++;
	int size_row = k - first_index;
	// if row not empty
	if (size_row > 0)
	  {
	    A.ReallocateRow(i,size_row);
	    k = first_index;
	    Sort(k, k+size_row-1, IndCol, Val);
	    for (int j = 0; j < size_row; j++)
	      {
		A.Index(i,j) = IndCol(k);
		A.Value(i,j) = Val(k);
		k++;
	      }
	  }
	else
	  A.ClearRow(i);
      }
  }
  
  
  //! Each row and column are scaled
  /*!
    We compute diag(S)*A*diag(S)
    where S = scale
  */
  template<class T1, class Allocator1, class T2, class Allocator2>
  void ScaleMatrix(const Vector<T2, Vect_Full, Allocator2>& scale,
		   Matrix<T1, General, ArrayRowSparse, Allocator1>& A)
  {
    int m = A.GetM();
    for (int i = 0; i < m; i++ )
      for (int j = 0; j < A.GetRowSize(i); j++ )
	A.Value(i,j) *= scale(i)*scale(A.Index(i,j));
    
  }
  
  
  //! Each row and column are scaled
  /*!
    We compute diag(S)*A*diag(S)
    where S = scale
  */
  template<class T1, class Allocator1, class T2, class Allocator2>
  void ScaleMatrix(const Vector<T2, Vect_Full, Allocator2>& scale,
		   Matrix<T1, Symmetric, ArrayRowSymSparse, Allocator1>& A)
  {
    int m = A.GetM();
    for (int i = 0; i < m; i++ )
      for (int j = 0; j < A.GetRowSize(i); j++ )
	A.Value(i,j) *= scale(i)*scale(A.Index(i,j));
    
  }
  
  
  //! Each row and column are scaled
  /*!
    We compute diag(S)*A*diag(S)
    where S = scale
  */
  template<class T1, class Allocator1, class T2, class Allocator2>
  void ScaleMatrix(const Vector<T2, Vect_Full, Allocator2>& scale,
		   Matrix<T1, Symmetric, ArrayRowSymComplexSparse,
		   Allocator1>& A)
  {
    int m = A.GetM();
    for (int i = 0; i < m; i++ )
      {
	for (int j = 0; j < A.GetRealRowSize(i); j++ )
	  A.ValueReal(i,j) *= scale(i)*scale(A.IndexReal(i,j));
	
	for (int j = 0; j < A.GetImagRowSize(i); j++ )
	  A.ValueImag(i,j) *= scale(i)*scale(A.IndexImag(i,j));
      }
    
  }
  
  
  //! Each row and column are scaled
  /*!
    We compute diag(S)*A*diag(S)
    where S = scale
  */
  template<class T1, class Allocator1, class T2, class Allocator2>
  void ScaleMatrix(const Vector<T2, Vect_Full, Allocator2>& scale,
		   Matrix<T1, General, ArrayRowComplexSparse, Allocator1>& A)
  {
    int m = A.GetM();
    for (int i = 0; i < m; i++ )
      {
	for (int j = 0; j < A.GetRealRowSize(i); j++ )
	  A.ValueReal(i,j) *= scale(i)*scale(A.IndexReal(i,j));
	
	for (int j = 0; j < A.GetImagRowSize(i); j++ )
	  A.ValueImag(i,j) *= scale(i)*scale(A.IndexImag(i,j));
      }
    
  }
  
  
  //! Each row is scaled
  /*!
    We compute diag(S)*A
    where S = scale
  */
  template<class T1, class Allocator1, class T2, class Allocator2>
  void ScaleLeftMatrix(const Vector<T2, Vect_Full, Allocator2>& scale,
		       Matrix<T1, General, ArrayRowSparse, Allocator1>& A)
  {
    int m = A.GetM();
    for (int i = 0; i < m; i++ )
      for (int j = 0; j < A.GetRowSize(i); j++ )
	A.Value(i,j) *= scale(i);
    
  }
  
  
  //! Each row is scaled
  /*!
    We compute diag(S)*A where S = scale.
    In order to keep symmetry, the operation
    is performed on upper part of the matrix,
    considering that lower part is affected by this
    operation
  */
  template<class T1, class Allocator1, class T2, class Allocator2>
  void ScaleLeftMatrix(const Vector<T2, Vect_Full, Allocator2>& scale,
		       Matrix<T1, Symmetric,
		       ArrayRowSymSparse, Allocator1>& A)
  {
    int m = A.GetM();
    for (int i = 0; i < m; i++ )
      for (int j = 0; j < A.GetRowSize(i); j++ )
	A.Value(i,j) *= scale(i);
    
  }
  
  
  //! Each row is scaled
  /*!
    We compute diag(S)*A where S = scale.
    In order to keep symmetry, the operation
    is performed on upper part of the matrix,
    considering that lower part is affected by this
    operation
  */
  template<class T1, class Allocator1, class T2, class Allocator2>
  void ScaleLeftMatrix(const Vector<T2, Vect_Full, Allocator2>& scale,
		       Matrix<T1, Symmetric,
		       ArrayRowSymComplexSparse, Allocator1>& A)
  {
    int m = A.GetM();
    for (int i = 0; i < m; i++ )
      {
	for (int j = 0; j < A.GetRealRowSize(i); j++ )
	  A.ValueReal(i,j) *= scale(i);
	
	for (int j = 0; j < A.GetImagRowSize(i); j++ )
	  A.ValueImag(i,j) *= scale(i);
      }
    
  }
  
  
  //! Each row is scaled
  /*!
    We compute diag(S)*A
    where S = scale
  */
  template<class T1, class Allocator1, class T2, class Allocator2>
  void ScaleLeftMatrix(const Vector<T2, Vect_Full, Allocator2>& scale,
		       Matrix<T1, General,
		       ArrayRowComplexSparse, Allocator1>& A)
  {
    int m = A.GetM();
    for (int i = 0; i < m; i++ )
      {
	for (int j = 0; j < A.GetRealRowSize(i); j++ )
	  A.ValueReal(i,j) *= scale(i);
	
	for (int j = 0; j < A.GetImagRowSize(i); j++ )
	  A.ValueImag(i,j) *= scale(i);
      }
    
  }
  
} // end namespace

#define SELDON_FILE_PERMUTATION_SCALING_MATRIX_CXX
#endif