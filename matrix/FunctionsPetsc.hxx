
#ifndef SELDON_FILE_FUNCTIONS_PETSC_HXX

namespace Seldon
{

  template <class T, class Allocator0, class Allocator1>
  void GetRow(const Matrix<T, General, PETScMPIAIJ, Allocator0>& M,
	      size_t i, Vector<T, PETScPar, Allocator1>& X);

  template <class T, class Allocator0, class Allocator1>
  void GetRow(const Matrix<T, General, PETScMPIDense, Allocator0>& M,
	      size_t i, Vector<T, PETScPar, Allocator1>& X);

  template <class T, class Allocator0,
            class Storage1, class Allocator1>
  void GetRow(const Matrix<T, General, PETScMPIDense, Allocator0>& M,
	      size_t i, Vector<T, Storage1, Allocator1>& X);

  template <class T0, class Prop0, class Allocator0,
	    class T1, class Allocator1>
  void GetCol(const Matrix<T0, Prop0, PETScSeqDense, Allocator0>& M,
        size_t j, Vector<T1, PETScSeq, Allocator1>& X);

  template <class T0, class Prop0, class Allocator0,
	    class T1, class Allocator1>
  void GetCol(const Matrix<T0, Prop0, PETScMPIDense, Allocator0>& M,
        size_t j, Vector<T1, PETScPar, Allocator1>& X);

  template <class T0, class Prop0, class Allocator0,
	    class T1, class Allocator1>
  void GetCol(const Matrix<T0, Prop0, PETScMPIAIJ, Allocator0>& M,
        size_t j, Vector<T1, PETScPar, Allocator1>& X);

  template <class T0, class Prop0, class Allocator0,
	    class T1, class Allocator1>
  void SetRow(const Vector<T1, PETScSeq, Allocator1>& X,
              size_t i, Matrix<T0, Prop0, PETScSeqDense, Allocator0>& M);

  template <class T0, class Prop0, class Allocator0,
	    class T1, class Allocator1>
  void SetRow(const Vector<T1, PETScPar, Allocator1>& X,
              size_t i, Matrix<T0, Prop0, PETScMPIDense, Allocator0>& M);

  template <class T0, class Prop0, class Allocator0,
	    class T1, class Allocator1>
  void SetRow(const Vector<T1, PETScPar, Allocator1>& X,
              size_t i, Matrix<T0, Prop0, PETScMPIAIJ, Allocator0>& M);

  template <class T0, class Prop0, class Allocator0,
            class T1, class Prop1, class Allocator1>
  void SetRow(const Vector<T1, Prop1, Allocator1>& X,
              size_t i, Matrix<T0, Prop0, PETScMPIDense, Allocator0>& M);

  template <class T0, class Prop0, class Allocator0,
	    class T1, class Allocator1>
  void SetCol(const Vector<T1, PETScSeq, Allocator1>& X,
              size_t j, Matrix<T0, Prop0, PETScSeqDense, Allocator0>& M);

  template <class T0, class Prop0, class Allocator0,
            class T1, class Prop1, class Allocator1>
  void SetCol(const Vector<T1, Prop1, Allocator1>& X,
              size_t j, Matrix<T0, Prop0, PETScMPIDense, Allocator0>& M);

  template <class T0, class Prop0, class Allocator0,
	    class T1, class Allocator1>
  void SetCol(const Vector<T1, PETScPar, Allocator1>& X,
              size_t j, Matrix<T0, Prop0, PETScMPIDense, Allocator0>& M);

  template <class T0, class Prop0, class Allocator0,
	    class T1, class Allocator1>
  void SetCol(const Vector<T1, PETScPar, Allocator1>& X,
              size_t j, Matrix<T0, Prop0, PETScMPIAIJ, Allocator0>& M);

  template <class T0, class Prop0, class Allocator0>
  void GetInverse(Matrix<T0, Prop0, PETScMPIDense, Allocator0>& M);

  template<class T, class Prop1, class Prop2, class Alloc1, class Alloc2>
  void Copy(const Matrix<T, Prop1, PETScMPIAIJ, Alloc1>& A,
		  Matrix<T, Prop2, PETScMPIAIJ, Alloc2>& B);

  template<class T, class Prop1, class Prop2, class Alloc1, class Alloc2>
  void Copy(const Matrix<T, Prop1, PETScMPIDense, Alloc1>& A,
		  Matrix<T, Prop2, PETScMPIDense, Alloc2>& B);

  template<class T, class Prop1, class Prop2, class Alloc1, class Alloc2>
  void Copy(const Matrix<T, Prop1, PETScMPIDense, Alloc1>& A,
		  Matrix<T, Prop2, PETScMPIAIJ, Alloc2>& B);

  template<class T, class Prop1, class Prop2, class Alloc1, class Alloc2>
  void Copy(const Matrix<T, Prop1, PETScMPIAIJ, Alloc1>& A,
		  Matrix<T, Prop2, PETScMPIDense, Alloc2>& B);

  template<class T, class Prop1, class Prop2,
           class Storage2, class Alloc1, class Alloc2>
  void Copy(const Matrix<T, Prop1, PETScMPIDense, Alloc1>& A,
		  Matrix<T, Prop2, Storage2, Alloc2>& B);


}


#define SELDON_FILE_FUNCTIONS_PETSC_HXX
#endif
