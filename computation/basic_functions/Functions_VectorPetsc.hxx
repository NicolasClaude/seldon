#ifndef SELDON_FILE_FUNCTIONS_VECTOR_PETSC_HXX
#define SELDON_FILE_FUNCTIONS_VECTOR_PETSC_HXX


namespace Seldon
{


  /////////
  // ADD //


  template <class T,
            class Allocator1,
            class Allocator2>
  void Add(const T alpha,
           const Vector<T, PETScSeq, Allocator1>& X,
           Vector<T, PETScSeq, Allocator2>& Y);

  template <class T,
            class Allocator1,
            class Allocator2>
  void Add(const T alpha,
           const Vector<T, PETScPar, Allocator1>& X,
           Vector<T, PETScPar, Allocator2>& Y);

  template <class T,
            class Allocator1,
            class Allocator2>
  void Add(const T alpha,
           const Vector<T, VectFull, Allocator1>& X,
           Vector<T, PETScPar, Allocator2>& Y);


  // ADD //
  /////////


} // namespace Seldon.


#endif
