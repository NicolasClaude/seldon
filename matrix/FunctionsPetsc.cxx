#ifndef SELDON_FILE_FUNCTIONS_PETSC_CXX

#include "FunctionsPetsc.hxx"


namespace Seldon
{


  //! Extracts a row from a PETSc matrix
  /*!
    \param M sparse matrix
    \param i row index
    \param X row extracted
    X = M(i, :)
  */
  template <class T, class Allocator0, class Allocator1>
  void GetRow(const Matrix<T, General, PETScMPIDense, Allocator0>& M,
	      size_t i, Vector<T, PETScPar, Allocator1>& X)
  {
#ifdef SELDON_CHECK_BOUNDS
    size_t m = M.GetM();
    if (i >= m)
      throw WrongIndex("GetRow()",
                       string("Index should be in [0, ") + to_str(m - 1)
                       + "], but is equal to " + to_str(i) + ".");
#endif
    X.Reallocate(M.GetN());
    int a, b;
    M.GetProcessorRowRange(a, b);

    if (i >= a && i < b)
      {
        int  ncols;
        const int * cols;
        const T * vals ;

        MatGetRow(M.GetPetscMatrix(), i, &ncols, &cols,  &vals);

        for (int j = 0; j < ncols ; j++ )
          X.SetBuffer(cols[j], vals[j]);

        MatRestoreRow(M.GetPetscMatrix(), i, &ncols, &cols,  &vals);
      }

    X.Flush();

  }

  //! Extracts a row from a PETSc matrix
  /*!
    \param M sparse matrix
    \param i row index
    \param X row extracted
    X = M(i, :)
  */
  template <class T, class Allocator0, class Allocator1>
  void GetRow(const Matrix<T, General, PETScMPIAIJ, Allocator0>& M,
	      size_t i, Vector<T, PETScPar, Allocator1>& X)
  {
#ifdef SELDON_CHECK_BOUNDS
    size_t m = M.GetM();
    if (i >= m)
      throw WrongIndex("GetRow()",
                       string("Index should be in [0, ") + to_str(m - 1)
                       + "], but is equal to " + to_str(i) + ".");
#endif
    X.Reallocate(M.GetN());
    int a, b;
    M.GetProcessorRowRange(a, b);

    if (i >= a && i < b)
      {
        int  ncols;
        const int * cols;
        const T * vals ;

        MatGetRow(M.GetPetscMatrix(), i, &ncols, &cols,  &vals);

        for (int j = 0; j < ncols ; j++ )
          X.SetBuffer(cols[j], vals[j]);

        MatRestoreRow(M.GetPetscMatrix(), i, &ncols, &cols,  &vals);
      }

    X.Flush();

  }


  //! Extracts a row from a PETSc matrix
  /*!
    \param M sparse matrix
    \param i row index
    \param X row extracted
    X = M(i, :)
  */
  template <class T, class Allocator0, class Storage1, class Allocator1>
  void GetRow(const Matrix<T, General, PETScMPIDense, Allocator0>& M,
	      size_t i, Vector<T, Storage1, Allocator1>& X)
  {
#ifdef SELDON_CHECK_BOUNDS
    size_t m = M.GetM();
    if (i >= m)
      throw WrongIndex("GetRow()",
                       string("Index should be in [0, ") + to_str(m - 1)
                       + "], but is equal to " + to_str(i) + ".");
#endif
    X.Reallocate(M.GetN());

    for (size_t j = 0; j < M.GetN(); j++)
      X(j) = M.GetOnAll(i, j);

  }




  //! Extracts a column from a matrix.
  /*!
    \param[in] M matrix.
    \param[in] j column index.
    \param[out] X extracted column.
  */
  template <class T0, class Prop0, class Allocator0,
            class T1, class Allocator1>
  void GetCol(const Matrix<T0, Prop0, PETScSeqDense, Allocator0>& M,
              size_t j, Vector<T1, PETScSeq, Allocator1>& X)
  {
    for (int i = 0; i < M.GetM(); i++)
      X.SetBuffer(i, M(i, j));
    X.Flush();
  }


  //! Extracts a column from a matrix.
  /*!
    \param[in] M matrix.
    \param[in] j column index.
    \param[out] X extracted column.
  */
  template <class T0, class Prop0, class Allocator0,
            class T1, class Allocator1>
  void GetCol(const Matrix<T0, Prop0, PETScMPIDense, Allocator0>& M,
              size_t j, Vector<T1, PETScPar, Allocator1>& X)
  {
    int a, b;
    M.GetProcessorRowRange(a, b);
    for (int i = a; i < b; i++)
      X.SetBuffer(i, M(i, j));
    X.Flush();
  }


  //! Extracts a column from a matrix.
  /*!
    \param[in] M matrix.
    \param[in] j column index.
    \param[out] X extracted column.
  */
  template <class T0, class Prop0, class Allocator0,
            class T1, class Allocator1>
  void GetCol(const Matrix<T0, Prop0, PETScMPIAIJ, Allocator0>& M,
              size_t j, Vector<T1, PETScPar, Allocator1>& X)
  {
    int a, b;
    M.GetProcessorRowRange(a, b);
    for (int i = a; i < b; i++)
      X.SetBuffer(i, M(i, j));
    X.Flush();
  }


  //! Sets a row of a matrix.
  /*!
    \param[in] X new row \a i of \a M.
    \param[in] i row index.
    \param[out] M matrix.
  */
  template <class T0, class Prop0, class Allocator0,
            class T1, class Allocator1>
  void SetRow(const Vector<T1, PETScSeq, Allocator1>& X,
              size_t i, Matrix<T0, Prop0, PETScSeqDense, Allocator0>& M)
  {
    for (size_t j = 0; j < M.GetN(); j++)
      M.SetBuffer(i, j, X(j));
    M.Flush();
  }


  //! Sets a row of a matrix.
  /*!
    \param[in] X new row \a i of \a M.
    \param[in] i row index.
    \param[out] M matrix.
  */
  template <class T0, class Prop0, class Allocator0,
            class T1, class Allocator1>
  void SetRow(const Vector<T1, PETScPar, Allocator1>& X,
              size_t i, Matrix<T0, Prop0, PETScMPIDense, Allocator0>& M)
  {
    for (size_t j = 0; j < M.GetN(); j++)
      M.SetBuffer(i, j, X.GetOnAll(j));
    M.Flush();
  }

  //! Sets a row of a matrix.
  /*!
    \param[in] X new row \a i of \a M.
    \param[in] i row index.
    \param[out] M matrix.
  */
  template <class T0, class Prop0, class Allocator0,
            class T1, class Allocator1>
  void SetRow(const Vector<T1, PETScPar, Allocator1>& X,
              size_t i, Matrix<T0, Prop0, PETScMPIAIJ, Allocator0>& M)
  {
    for (size_t j = 0; j < M.GetN(); j++)
      M.SetBuffer(i, j, X.GetOnAll(j));
    M.Flush();
  }

  //! Sets a row of a matrix.
  /*!
    \param[in] X new row \a i of \a M.
    \param[in] i row index.
    \param[out] M matrix.
  */
  template <class T0, class Prop0, class Allocator0,
            class T1, class Prop1, class Allocator1>
  void SetRow(const Vector<T1, Prop1, Allocator1>& X,
              size_t i, Matrix<T0, Prop0, PETScMPIAIJ, Allocator0>& M)
  {
    for (size_t j = 0; j < M.GetN(); j++)
      M.SetBuffer(i, j, X(j));
    M.Flush();
  }


 //! Sets a column of a matrix.
  /*!
    \param[in] X new column \a j of \a M.
    \param[in] j column index.
    \param[out] M matrix.
  */
  template <class T0, class Prop0, class Allocator0,
            class T1, class Allocator1>
  void SetCol(const Vector<T1, PETScSeq, Allocator1>& X,
              size_t j, Matrix<T0, Prop0, PETScSeqDense, Allocator0>& M)
  {
    int a, b;
    M.GetProcessorRowRange(a, b);
    for (int i = a; i < b; i++)
      M.SetBuffer(i, j, X(i));
    M.Flush();
  }


  //! Sets a column of a matrix.
  /*!
    \param[in] X new column \a j of \a M.
    \param[in] j column index.
    \param[out] M matrix.
  */
  template <class T0, class Prop0, class Allocator0,
            class T1, class Allocator1>
  void SetCol(const Vector<T1, PETScPar, Allocator1>& X,
              size_t j, Matrix<T0, Prop0, PETScMPIDense, Allocator0>& M)
  {
    int a, b;
    M.GetProcessorRowRange(a, b);
    for (int i = a; i < b; i++)
      M.SetBuffer(i, j, X(i));
    M.Flush();
  }


  //! Sets a column of a matrix.
  /*!
    \param[in] X new column \a j of \a M.
    \param[in] j column index.
    \param[out] M matrix.
  */
  template <class T0, class Prop0, class Allocator0,
            class T1, class Allocator1>
  void SetCol(const Vector<T1, PETScPar, Allocator1>& X,
              size_t j, Matrix<T0, Prop0, PETScMPIAIJ, Allocator0>& M)
  {
    int a, b;
    M.GetProcessorRowRange(a, b);
    for (int i = a; i < b; i++)
      M.SetBuffer(i, j, X(i));
    M.Flush();
  }


  //! Sets a column of a matrix.
  /*!
    \param[in] X new column \a j of \a M.
    \param[in] j column index.
    \param[out] M matrix.
  */
  template <class T0, class Prop0, class Allocator0,
            class T1, class Prop1, class Allocator1>
  void SetCol(const Vector<T1, Prop1, Allocator1>& X,
              size_t j, Matrix<T0, Prop0, PETScMPIDense, Allocator0>& M)
  {
    int a, b;
    M.GetProcessorRowRange(a, b);
    for (int i = a; i < b; i++)
      M.SetBuffer(i, j, X(i));
    M.Flush();
  }


  template <class T0, class Prop0, class Allocator0>
  void GetInverse(Matrix<T0, Prop0, PETScMPIDense, Allocator0>& M)
  {
    Matrix<T0> temp;
    Copy(M, temp);
    GetInverse(temp);
    Copy(temp, M);
  }


  template<class T, class Prop1, class Prop2, class Alloc1, class Alloc2>
  void Copy(const Matrix<T, Prop1, PETScMPIAIJ, Alloc1>& A,
		  Matrix<T, Prop2, PETScMPIAIJ, Alloc2>& B)
  {
    B.Copy(A);
  }


  template<class T, class Prop1, class Prop2, class Alloc1, class Alloc2>
  void Copy(const Matrix<T, Prop1, PETScMPIAIJ, Alloc1>& A,
		  Matrix<T, Prop2, PETScMPIDense, Alloc2>& B)
  {
    int a,b;
    B.GetProcessorRowRange(a, b);
    for (size_t i = a; i < b; i++)
      for (size_t j = 0; j < A.GetN(); j++)
        B.SetBuffer(i, j, A(i,j));
    B.Flush();
  }


  template<class T, class Prop1, class Prop2, class Alloc1, class Alloc2>
  void Copy(const Matrix<T, Prop1, RowMajor, Alloc1>& A,
		  Matrix<T, Prop2, PETScMPIDense, Alloc2>& B)
  {
    B.Reallocate(A.GetM(), A.GetN());
    int a,b;
    B.GetProcessorRowRange(a, b);
    for (size_t i = a; i < b; i++)
      for (size_t j = 0; j < A.GetN(); j++)
        B.SetBuffer(i, j, A(i,j));
    B.Flush();
  }


  template<class T, class Prop1, class Prop2, class Alloc1, class Alloc2>
  void Copy(const Matrix<T, Prop1, PETScMPIDense, Alloc1>& A,
		  Matrix<T, Prop2, PETScMPIAIJ, Alloc2>& B)
  {
    int a,b;
    B.GetProcessorRowRange(a, b);
    for (size_t i = a; i < b; i++)
      for (size_t j = 0; j < A.GetN(); j++)
        B.SetBuffer(i, j, A(i,j));
    B.Flush();
  }


  template<class T, class Prop1, class Prop2, class Alloc1, class Alloc2>
  void Copy(const Matrix<T, Prop1, PETScMPIDense, Alloc1>& A,
		  Matrix<T, Prop2, PETScMPIDense, Alloc2>& B)
  {
    B.Copy(A);
  }

  template<class T, class Prop1, class Prop2,
           class Storage2, class Alloc1, class Alloc2>
  void Copy(const Matrix<T, Prop1, PETScMPIDense, Alloc1>& A,
		  Matrix<T, Prop2, Storage2, Alloc2>& B)
  {
    int low, high, m, n;
    m = A.GetM();
    n = A.GetN();
    B.Reallocate(m, n);
    A.GetProcessorRowRange(low, high);
    int nrows = high - low;
    int rows[nrows];
    int cols[n];
    T values[m * n] = {0};

    for (int i = 0; i < nrows; i++)
      rows[i] = low + i;
    for (int i = 0; i < n; i++)
      cols[i] = i;

    MatGetValues(A.GetPetscMatrix(), nrows, rows, n,
                 cols, &values[low * n]);

    typename Alloc2::pointer reduced = Alloc2::allocate(n * m);


    MPI_Allreduce(values, reduced, m * n ,
                  MPI_DOUBLE, MPI_SUM, A.GetCommunicator());


    B.SetData(m, n , reduced);
  }


} // namespace Seldon


#define SELDON_FILE_FUNCTIONS_PETSC_CXX
#endif
