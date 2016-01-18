// Copyright (C) 2001-2011 Vivien Mallet
// Copyright (C) 2003-2009 Marc Duruflé
// Copyright (C) 2010 INRIA
// Author(s): Nicolas Claude
//
// This file is part of the linear-algebra library Seldon,
// http://seldon.sourceforge.net/.
//
// Seldon is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 2.1 of the License, or (at your option)
// any later version.
//
// Seldon is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
// more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Seldon. If not, see http://www.gnu.org/licenses/.


#ifndef SELDON_FILE_FUNCTIONS_MATRIX_PETSC_HXX
#define SELDON_FILE_FUNCTIONS_MATRIX_PETSC_HXX


namespace Seldon
{


  /////////
  // ADD //


  template <class T,
	    class Prop1, class Allocator1,
            class Prop2, class Allocator2>
  void AddMatrix(const T& alpha,
                 const Matrix<T, Prop1, PETScMPIAIJ, Allocator1>& A,
                 Matrix<T, Prop2, PETScMPIAIJ, Allocator2>& B);

  template <class T,
	    class Prop1, class Allocator1,
            class Prop2, class Allocator2>
  void AddMatrix(const T& alpha,
                 const Matrix<T, Prop1, PETScMPIDense, Allocator1>& A,
                 Matrix<T, Prop2, PETScMPIDense, Allocator2>& B);


  // ADD //
  /////////


  /////////
  // MLT //


  template <class T,
	    class Prop1, class Allocator1>
  void MltScalar(const T alpha,
                 Matrix<T, Prop1, PETScMPIAIJ, Allocator1>& A);

  template <class T,
	    class Prop1, class Allocator1>
  void MltScalar(const T alpha,
                 Matrix<T, Prop1, PETScMPIDense, Allocator1>& A);



  // MLT //
  /////////


  ////////////
  // MLTADD //


  template <class T0,
            class T1, class Prop1, class Allocator1,
            class T2, class Allocator2,
            class T3,
            class T4, class Prop4, class Allocator4>
  void MltAdd(const T0 alpha,
              const Matrix<T1, Prop1, PETScMPIDense, Allocator1>& A,
              const Matrix<T2, General, RowMajor, Allocator2>& B,
              const T3 beta,
              Matrix<T4, Prop4, PETScMPIDense, Allocator4>& C);

  template <class T,
            class Prop1, class Allocator1,
            class Allocator2,
            class Prop4, class Allocator4>
  void MltAdd(const T alpha,
              const Matrix<T, Prop1, PETScMPIAIJ, Allocator1>& A,
              const Matrix<T, General, RowMajor, Allocator2>& B,
              const T beta,
              Matrix<T, Prop4, PETScMPIAIJ, Allocator4>& C);

  template <class T,
            class Prop1, class Allocator1,
            class Allocator2,
            class Prop4, class Allocator4>
  void MltAdd(const T alpha,
              const Matrix<T, Prop1, PETScMPIAIJ, Allocator1>& A,
              const Matrix<T, General, RowMajor, Allocator2>& B,
              const T beta,
              Matrix<T, Prop4, PETScMPIAIJ, Allocator4>& C);

  template <class T,
            class Prop1, class Allocator1,
            class Prop2, class Allocator2,
            class Prop4, class Allocator4>
  void MltAdd(const T alpha,
              const Matrix<T, Prop1, PETScMPIAIJ, Allocator1>& A,
              const Matrix<T, Prop2, PETScMPIAIJ, Allocator2>& B,
              const T beta,
              Matrix<T, Prop4, PETScMPIAIJ, Allocator4>& C);

  template <class T,
            class Prop1, class Allocator1,
            class Prop2, class Allocator2,
            class Prop4, class Allocator4>
  void MltAdd(const T alpha,
              const SeldonTranspose& TransA,
              const Matrix<T, Prop1, PETScMPIAIJ, Allocator1>& A,
              const SeldonTranspose& TransB,
              const Matrix<T, Prop2, PETScMPIAIJ, Allocator2>& B,
              const T beta,
              Matrix<T, Prop4, PETScMPIAIJ, Allocator4>& C);

  template <class T,
            class Prop1, class Allocator1,
            class Prop2, class Allocator2,
            class Prop4, class Allocator4>
  void MltAddMatrix(const T alpha,
              const Matrix<T, Prop1, PETScMPIDense, Allocator1>& A,
              const Matrix<T, Prop2, PETScMPIDense, Allocator2>& B,
              const T beta,
              Matrix<T, Prop4, PETScMPIDense, Allocator4>& C);

  template <class T,
            class Prop1, class Allocator1,
            class Prop2, class Allocator2,
            class Prop4, class Allocator4>
  void MltAddMatrix(const T alpha,
              const SeldonTranspose& TransA,
              const Matrix<T, Prop1, PETScMPIDense, Allocator1>& A,
              const SeldonTranspose& TransB,
              const Matrix<T, Prop2, PETScMPIDense, Allocator2>& B,
              const T beta,
              Matrix<T, Prop4, PETScMPIDense, Allocator4>& C);

  template <class T,
	    class Prop1, class Allocator1,
	    class Prop2, class Storage2, class Allocator2,
	    class Prop4, class Storage4, class Allocator4>
  void MltAdd(const T alpha,
	      const SeldonTranspose& TransA,
	      const Matrix<T, Prop1, PETScMPIAIJ, Allocator1>& A,
	      const SeldonTranspose& TransB,
	      const Matrix<T, Prop2, Storage2, Allocator2>& B,
	      const T beta,
	      Matrix<T, Prop4, Storage4, Allocator4>& C);
  ///////////////
  // TRANSPOSE //


  template<class T, class Prop, class Allocator>
  void Transpose(Matrix<T, Prop, PETScMPIAIJ, Allocator>& A);


  template<class T, class Prop, class Allocator>
  void Transpose(Matrix<T, Prop, PETScMPIDense, Allocator>& A);


  // TRANSPOSE //
  ///////////////


  ///////////
  // GETLU //


  template <class T, class Prop0, class Allocator0>
  void GetLU(Matrix<T, Prop0, PETScMPIDense, Allocator0>& A);


  // GETLU //
  ///////////


} // namespace seldon


#endif
