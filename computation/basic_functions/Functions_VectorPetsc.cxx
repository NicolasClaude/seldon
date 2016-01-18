
#ifndef SELDON_FILE_FUNCTIONS_VECTOR_PETSC_CXX
#define SELDON_FILE_FUNCTIONS_VECTOR_PETSC_CXX

#include "Functions_VectorPetsc.hxx"



namespace Seldon
{


  /////////
  // ADD //



  template <class T,
            class Allocator1,
            class Allocator2>
  void Add(const T alpha,
           const Vector<T, PETScSeq, Allocator1>& X,
           Vector<T, PETScSeq, Allocator2>& Y)
  {
    if (alpha != T(0))
      {
#ifdef SELDON_CHECK_DIMENSIONS
        CheckDim(X, Y, "Add(alpha, X, Y)");
#endif

        VecAXPY(Y.GetPetscVector(), alpha, X.GetPetscVector());
      }
  }


  template <class T,
            class Allocator1,
            class Allocator2>
  void Add(const T alpha,
           const Vector<T, PETScPar, Allocator1>& X,
           Vector<T, PETScPar, Allocator2>& Y)
  {
    if (alpha != T(0))
      {
#ifdef SELDON_CHECK_DIMENSIONS
        CheckDim(X, Y, "Add(alpha, X, Y)");
#endif

        VecAXPY(Y.GetPetscVector(), alpha, X.GetPetscVector());
      }
  }


  template <class T,
            class Allocator1,
            class Allocator2>
  void Add(const T alpha,
           const Vector<T, VectFull, Allocator1>& X,
           Vector<T, PETScPar, Allocator2>& Y)
  {
    if (alpha != T(0))
      {
#ifdef SELDON_CHECK_DIMENSIONS
        CheckDim(X, Y, "Add(alpha, X, Y)");
#endif

        int a,b;
        Y.GetProcessorRange(a, b);
        for (int i = a; i < b; i++)
          Y.SetBuffer(i, Y(i) + alpha * X(i));
        Y.Flush();
      }
  }


  // ADD //
  /////////


  //////////
  // COPY //


  template<class T, class Alloc1, class Alloc2>
  void Copy(const Vector<T, PETScPar, Alloc1>& A,
            Vector<T, VectFull, Alloc2>& B)
  {
    B.Reallocate(A.GetSize());
    for (int i = 0; i < A.GetSize(); i++)
      B(i) = A.GetOnAll(i);
  }


  template<class T, class Alloc1, class Alloc2>
  void Copy(const Vector<T, VectFull, Alloc1>& A,
            Vector<T, PETScPar, Alloc2>& B)
  {
    B.Reallocate(A.GetSize());
    int low, high;
    B.GetProcessorRange(low, high);
    for (int i = low; i < high; i++)
      B.SetBuffer(i, A(i));
    B.Flush();
  }


  // COPY //
  //////////


} // namespace Seldon.


#endif
