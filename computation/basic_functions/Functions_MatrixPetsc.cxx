#ifndef SELDON_FILE_FUNCTIONS_MATRIX_PETSC_CXX
#define SELDON_FILE_FUNCTIONS_MATRIX_PETSC_CXX


#include "Functions_MatrixPetsc.hxx"


namespace Seldon
{


  /////////
  // ADD //


  template <class T,
	    class Prop1, class Allocator1,
            class Prop2, class Allocator2>
  void AddMatrix(const T& alpha,
                 const Matrix<T, Prop1, PETScMPIAIJ, Allocator1>& A,
                 Matrix<T, Prop2, PETScMPIAIJ, Allocator2>& B)
  {
    MatAXPY(B.GetPetscMatrix(), alpha,
            A.GetPetscMatrix(), DIFFERENT_NONZERO_PATTERN);
  }


  template <class T,
	    class Prop1, class Allocator1,
            class Prop2, class Allocator2>
  void AddMatrix(const T& alpha,
                 const Matrix<T, Prop1, PETScMPIDense, Allocator1>& A,
                 Matrix<T, Prop2, PETScMPIDense, Allocator2>& B)
  {
    MatAXPY(B.GetPetscMatrix(), alpha,
            A.GetPetscMatrix(), DIFFERENT_NONZERO_PATTERN);
  }


  // ADD //
  /////////


  /////////
  // MLT //


  //! Multiplies a matrix by a scalar.
  /*!
    \param[in] alpha scalar.
    \param[in,out] M matrix to be multiplied.
  */
  template <class T,
	    class Prop1, class Allocator1>
  void MltScalar(const T alpha, Matrix<T, Prop1, PETScMPIAIJ, Allocator1>& A)
  {
    MatScale(A.GetPetscMatrix(), alpha);
  }


  //! Multiplies a matrix by a scalar.
  /*!
    \param[in] alpha scalar.
    \param[in,out] M matrix to be multiplied.
  */
  template <class T,
	    class Prop1, class Allocator1>
  void MltScalar(const T alpha, Matrix<T, Prop1, PETScMPIDense, Allocator1>& A)
  {
    MatScale(A.GetPetscMatrix(), alpha);
  }


  // MLT //
  /////////


  ////////////
  // MLTADD //


//! Multiplies two matrices, and adds the result to a third matrix.
  /*! It performs the operation \f$ C = \alpha A B + \beta C \f$ where \f$
    \alpha \f$ and \f$ \beta \f$ are scalars, and \f$ A \f$, \f$ B \f$ and \f$
    C \f$ are matrices.
    \param[in] alpha scalar.
    \param[in] A matrix.
    \param[in] B matrix.
    \param[in] beta scalar.
    \param[in,out] C matrix, result of the product of \a A with \a B, times \a
    alpha, plus \a beta times \a C.
  */
  template <class T,
            class Prop1, class Allocator1,
            class Allocator2,
            class Prop4, class Allocator4>
  void MltAdd(const T alpha,
              const Matrix<T, Prop1, PETScMPIDense, Allocator1>& A,
              const Matrix<T, General, RowMajor, Allocator2>& B,
              const T beta,
              Matrix<T, Prop4, PETScMPIDense, Allocator4>& C)
  {
    int na = A.GetN();
    int mc = C.GetM();
    int nc = C.GetN();

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, B, C, "MltAdd(alpha, A, B, beta, C)");
#endif
    T *local_a;
    MatGetArray(A.GetPetscMatrix(), &local_a);
    int nlocal_A;
    int mlocal_A;
    MatGetLocalSize(A.GetPetscMatrix(), &mlocal_A, &nlocal_A);
    Matrix<T, Prop1, ColMajor, Allocator1> local_A;
    local_A.SetData(mlocal_A, na, local_a);

    T *local_c;
    MatGetArray(C.GetPetscMatrix(), &local_c);
    int nlocal_C;
    int mlocal_C;
    MatGetLocalSize(C.GetPetscMatrix(), &mlocal_C, &nlocal_C);
    Matrix<T, Prop4, ColMajor, Allocator4> local_C;
    local_C.SetData(mlocal_C, nc, local_c);

    MltAdd(alpha, local_A, B, beta, local_C);

    local_A.Nullify();
    MatRestoreArray(A.GetPetscMatrix(), &local_a);
    A.Flush();

    local_C.Nullify();
    MatRestoreArray(C.GetPetscMatrix(), &local_c);
    C.Flush();
  }


  //! Multiplies two matrices, and adds the result to a third matrix.
  /*! It performs the operation \f$ C = \alpha A B + \beta C \f$ where \f$
    \alpha \f$ and \f$ \beta \f$ are scalars, and \f$ A \f$, \f$ B \f$ and \f$
    C \f$ are matrices.
    \param[in] alpha scalar.
    \param[in] A matrix.
    \param[in] B matrix.
    \param[in] beta scalar.
    \param[in,out] C matrix, result of the product of \a A with \a B, times \a
    alpha, plus \a beta times \a C.
  */
  template <class T,
            class Prop1, class Allocator1,
            class Allocator2,
            class Prop4, class Allocator4>
  void MltAdd(const T alpha,
              const Matrix<T, Prop1, PETScMPIAIJ, Allocator1>& A,
              const Matrix<T, General, RowMajor, Allocator2>& B,
              const T beta,
              Matrix<T, Prop4, PETScMPIAIJ, Allocator4>& C)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, B, C, "MltAdd(alpha, A, B, beta, C)");
#endif

    Matrix<T> tmp;
    Copy(B, tmp);
    MltAdd(alpha, A, tmp, beta, C);
  }


  //! Multiplies two matrices, and adds the result to a third matrix.
  /*! It performs the operation \f$ C = \alpha A B + \beta C \f$ where \f$
    \alpha \f$ and \f$ \beta \f$ are scalars, and \f$ A \f$, \f$ B \f$ and \f$
    C \f$ are matrices.
    \param[in] alpha scalar.
    \param[in] A matrix.
    \param[in] B matrix.
    \param[in] beta scalar.
    \param[in,out] C matrix, result of the product of \a A with \a B, times \a
    alpha, plus \a beta times \a C.
  */
  template <class T,
            class Prop1, class Allocator1,
            class Allocator2,
            class Allocator4>
  void MltAdd(const T alpha,
              const Matrix<T, Prop1, PETScMPIAIJ, Allocator1>& A,
              const Matrix<T, General, RowMajor, Allocator2>& B,
              const T beta,
              Matrix<T, General, RowMajor, Allocator4>& C)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, B, C, "MltAdd(alpha, A, B, beta, C)");
#endif

    Matrix<T> tmp;
    Copy(A, tmp);
    MltAdd(alpha, tmp, B, beta, C);
  }


  //! Multiplies two matrices, and adds the result to a third matrix.
  /*! It performs the operation \f$ C = \alpha A B + \beta C \f$ where \f$
    \alpha \f$ and \f$ \beta \f$ are scalars, and \f$ A \f$, \f$ B \f$ and \f$
    C \f$ are matrices.
    \param[in] alpha scalar.
    \param[in] A matrix.
    \param[in] B matrix.
    \param[in] beta scalar.
    \param[in,out] C matrix, result of the product of \a A with \a B, times \a
    alpha, plus \a beta times \a C.
  */
  template <class T,
            class Prop1, class Allocator1,
            class Prop2, class Allocator2,
            class Prop4, class Allocator4>
  void MltAdd(const T alpha,
              const Matrix<T, Prop1, PETScMPIAIJ, Allocator1>& A,
              const Matrix<T, Prop2, PETScMPIAIJ, Allocator2>& B,
              const T beta,
              Matrix<T, Prop4, PETScMPIAIJ, Allocator4>& C)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, B, C, "MltAdd(alpha, A, B, beta, C)");
#endif


    MatScale(C.GetPetscMatrix(), beta);

    Mat tmp;
    MatMatMultSymbolic(A.GetPetscMatrix(),B.GetPetscMatrix(),
                       PETSC_DEFAULT, &tmp);

    MatMatMultNumeric(A.GetPetscMatrix(),B.GetPetscMatrix(),
                      tmp);
    // C = C + alpha * temp

    MatAXPY(C.GetPetscMatrix(), alpha,
            tmp, DIFFERENT_NONZERO_PATTERN);
    MatDestroy(&tmp);

  }


  template <class T,
            class Prop1, class Allocator1,
            class Prop2, class Allocator2,
            class Prop4, class Allocator4>
  void MltAddMatrix(const T alpha,
              const Matrix<T, Prop1, PETScMPIDense, Allocator1>& A,
              const Matrix<T, Prop2, PETScMPIDense, Allocator2>& B,
              const T beta,
              Matrix<T, Prop4, PETScMPIDense, Allocator4>& C)
  {
#ifdef SELDON_CHECK_DIMENSIONS
        CheckDim(A, B, C, "MltAdd(alpha,  A,  B, beta, C)");
#endif
    if (alpha == 0)
      {
        MatScale(C.GetPetscMatrix(), beta);
        return;
      }
    if (beta != T(1))
      MatScale(C.GetPetscMatrix(), beta);

    Matrix<T> copyB;
    Copy(B, copyB);

    int high,low;
    A.GetProcessorRowRange(low, high);
    int nrows = high - low;
    T* dataA = new T[nrows * A.GetN()];
    A.GetLocalData(dataA);
    T* dataB  = copyB.GetData();

    for (int i = 0; i < nrows; i++)
      for (int j = 0; j < A.GetN(); j++)
        {
          T value = 0;
          for (int k = 0; k < A.GetN(); k++)
            value += dataA[i * A.GetN() + k] * dataB[k * B.GetN() + j];
          C.SetBuffer(i + low, j, alpha * value, ADD_VALUES);
        }
    C.Flush();
    delete[] dataA;
  }


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
              Matrix<T, Prop4, PETScMPIDense, Allocator4>& C)
  {
    if (!TransA.NoTrans())
      throw WrongArgument("MltAdd(const T0 alpha, SeldonTranspose TransA, "
                          "const Matrix<T1, Prop1, Storage1, Allocator1>& A"
                          "SeldonTranspose TransB,"
                          "const Matrix<T2, Prop2, Storage2, Allocator2>& B,"
                          "const T3 beta,"
                          "Matrix<T4, Prop4, Storage4, Allocator4>& C)",
                          "'TransA' must be equal to 'SeldonNoTrans'.");
    if (!TransB.NoTrans() && !TransB.Trans())
      throw WrongArgument("MltAdd(const T0 alpha, SeldonTranspose TransA, "
                          "const Matrix<T1, Prop1, Storage1, Allocator1>& A"
                          "SeldonTranspose TransB,"
                          "const Matrix<T2, Prop2, Storage2, Allocator2>& B,"
                          "const T3 beta,"
                          "Matrix<T4, Prop4, Storage4, Allocator4>& C)",
                          "'TransB' must be equal to 'SeldonNoTrans' or "
                          "'SeldonTrans'.");

    if (!TransB.Trans())
      MltAdd(alpha, A, B, beta, C);
    else
      {
#ifdef SELDON_CHECK_DIMENSIONS
        CheckDim(SeldonNoTrans, A, SeldonTrans, B, C,
                 "MltAdd(alpha, TransA, A, TransB, B, beta, C)");
#endif

        if (alpha == 0)
          {
            MatScale(C.GetPetscMatrix(), beta);
            return;
          }
        if (beta != T(1))
          MatScale(C.GetPetscMatrix(), beta);

        Matrix<T> copyB;
        Copy(B, copyB);

        int high,low;
        A.GetProcessorRowRange(low, high);
        int nrows = high - low;
        T* dataA = new T[nrows * A.GetN()];
        A.GetLocalData(dataA);
        T* dataB = copyB.GetData();

        for (int i = 0; i < nrows; i++)
          {
            for (int j = 0; j < A.GetN(); j++)
              {
                T value = 0;
                for (size_t k = 0; k < A.GetN(); k++)
                  value += dataA[i * A.GetN() + k] * dataB[j * B.GetN() + k];

                C.SetBuffer(i + low, j, alpha * value, ADD_VALUES);
              }
          }
        C.Flush();
        delete[] dataA;
      }

  }


  //! Multiplies two matrices, and adds the result to a third matrix.
  /*! It performs the operation \f$ C = \alpha A B + \beta C \f$
    or \f$C = \alpha A B^T + \beta C \f$, where \f$
    \alpha \f$ and \f$ \beta \f$ are scalars, and \f$ A \f$, \f$ B \f$ and \f$
    C \f$ are matrices.
    \param[in] alpha scalar.
    \param[in] TransA status of A: it must be SeldonNoTrans. This argument
    is required for consistency with the interface for full matrices.
    \param[in] A matrix.
    \param[in] TransB status of B: SeldonNoTrans or SeldonTrans.
    \param[in] B matrix.
    \param[in] beta scalar.
    \param[in,out] C matrix, result of the product of \a A with \a B or \a
    B^T, times \a alpha, plus \a beta times \a C.
  */
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
	      Matrix<T, Prop4, PETScMPIAIJ, Allocator4>& C)
  {
    if (!TransA.NoTrans())
      throw WrongArgument("MltAdd(const T0 alpha, SeldonTranspose TransA, "
                          "const Matrix<T1, Prop1, Storage1, Allocator1>& A"
                          "SeldonTranspose TransB,"
                          "const Matrix<T2, Prop2, Storage2, Allocator2>& B,"
                          "const T3 beta,"
                          "Matrix<T4, Prop4, Storage4, Allocator4>& C)",
                          "'TransA' must be equal to 'SeldonNoTrans'.");
    if (!TransB.NoTrans() && !TransB.Trans())
      throw WrongArgument("MltAdd(const T0 alpha, SeldonTranspose TransA, "
                          "const Matrix<T1, Prop1, Storage1, Allocator1>& A"
                          "SeldonTranspose TransB,"
                          "const Matrix<T2, Prop2, Storage2, Allocator2>& B,"
                          "const T3 beta,"
                          "Matrix<T4, Prop4, Storage4, Allocator4>& C)",
                          "'TransB' must be equal to 'SeldonNoTrans' or "
                          "'SeldonTrans'.");

    if (!TransB.Trans())
      MltAdd(alpha, A, B, beta, C);
    else
      {

#ifdef SELDON_CHECK_DIMENSIONS
        CheckDim(SeldonNoTrans, A, SeldonTrans, B, C,
                 "MltAdd(alpha, TransA, A, TransB, B, beta, C)");
#endif
        int m,n;
        Mat* sub_B;
        IS rowIS,colIS;
        ISCreateStride(MPI_COMM_SELF, B.GetM(), 0, 1, &rowIS);
        ISCreateStride(MPI_COMM_SELF, B.GetN(), 0, 1, &colIS);

        MatGetOwnershipRange(A.GetPetscMatrix(),&m,&n);
        MatGetSubMatrices(B.GetPetscMatrix(), 1,
                          const_cast<IS*>(&rowIS),
                          const_cast<IS*>(&colIS),
                          MAT_INITIAL_MATRIX, &sub_B);


        for (unsigned int i = m; i < n; i++)
          {
            int non_zero_A, non_zero_B;
            const int* col_A;
            const int* col_B;
            const double* value_A;
            const double* value_B;
            MatGetRow(A.GetPetscMatrix(), i, &non_zero_A, &col_A, &value_A);
            for (unsigned int j = 0; j < A.GetN(); j++)
              {
                T value = 0;
                MatGetRow(sub_B[0], j, &non_zero_B, &col_B, &value_B);
                int
                  iterator_A = 0,
                  iterator_B = 0;
                while (iterator_A < non_zero_A && iterator_B < non_zero_B)
                  {
                    if (col_A[iterator_A] == col_B[iterator_B])
                      {
                        value += value_A[iterator_A] * value_B[iterator_B];
                        ++iterator_A;
                        ++iterator_B;
                      }
                    else if (col_A[iterator_A] < col_B[iterator_B])
                      ++iterator_A;
                    else
                      ++iterator_B;

                  }//while
                MatRestoreRow(sub_B[0], i, &non_zero_B, &col_B, &value_B);
                value *= alpha;
                if (beta != 0)
                  value += C(i,j) * beta;
                C.SetBuffer(i, j, value);
              } // j
            MatRestoreRow(A.GetPetscMatrix(), i,
                          &non_zero_A, &col_A, &value_A);
          } // boucle i

        C.Flush();
        MatDestroyMatrices(1, &sub_B);
        ISDestroy(&rowIS);
        ISDestroy(&colIS);
      }
  }


 //! Multiplies two matrices, and adds the result to a third matrix.
  /*! It performs the operation \f$ C = \alpha A B + \beta C \f$
    or \f$C = \alpha A B^T + \beta C \f$, where \f$
    \alpha \f$ and \f$ \beta \f$ are scalars, and \f$ A \f$, \f$ B \f$ and \f$
    C \f$ are matrices.
    \param[in] alpha scalar.
    \param[in] TransA status of A: it must be SeldonNoTrans. This argument
    is required for consistency with the interface for full matrices.
    \param[in] A matrix.
    \param[in] TransB status of B: SeldonNoTrans or SeldonTrans.
    \param[in] B matrix.
    \param[in] beta scalar.
    \param[in,out] C matrix, result of the product of \a A with \a B or \a
    B^T, times \a alpha, plus \a beta times \a C.
  */
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
	      Matrix<T, Prop4, Storage4, Allocator4>& C)
  {
    if (!TransA.NoTrans())
      throw WrongArgument("MltAdd(const T0 alpha, SeldonTranspose TransA, "
                          "const Matrix<T1, Prop1, Storage1, Allocator1>& A"
                          "SeldonTranspose TransB,"
                          "const Matrix<T2, Prop2, Storage2, Allocator2>& B,"
                          "const T3 beta,"
                          "Matrix<T4, Prop4, Storage4, Allocator4>& C)",
                          "'TransA' must be equal to 'SeldonNoTrans'.");
    if (!TransB.NoTrans() && !TransB.Trans())
      throw WrongArgument("MltAdd(const T0 alpha, SeldonTranspose TransA, "
                          "const Matrix<T1, Prop1, Storage1, Allocator1>& A"
                          "SeldonTranspose TransB,"
                          "const Matrix<T2, Prop2, Storage2, Allocator2>& B,"
                          "const T3 beta,"
                          "Matrix<T4, Prop4, Storage4, Allocator4>& C)",
                          "'TransB' must be equal to 'SeldonNoTrans' or "
                          "'SeldonTrans'.");

    if (!TransB.Trans())
      MltAdd(alpha, A, B, beta, C);
    else
      {

#ifdef SELDON_CHECK_DIMENSIONS
        CheckDim(SeldonNoTrans, A, SeldonTrans, B, C,
                 "MltAdd(alpha, TransA, A, TransB, B, beta, C)");
#endif


        Matrix<T> tmp;
        Copy(A, tmp);
        MltAdd(alpha, TransA, tmp, TransB, B, beta, C);
      }
  }




//! Multiplies two matrices, and adds the result to a third matrix.
  /*! It performs the operation \f$ C = \alpha A B + \beta C \f$
    or \f$C = \alpha A B^T + \beta C \f$, where \f$
    \alpha \f$ and \f$ \beta \f$ are scalars, and \f$ A \f$, \f$ B \f$ and \f$
    C \f$ are matrices.
    \param[in] alpha scalar.
    \param[in] TransA status of A: it must be SeldonNoTrans. This argument
    is required for consistency with the interface for full matrices.
    \param[in] A matrix.
    \param[in] TransB status of B: SeldonNoTrans or SeldonTrans.
    \param[in] B matrix.
    \param[in] beta scalar.
    \param[in,out] C matrix, result of the product of \a A with \a B or \a
    B^T, times \a alpha, plus \a beta times \a C.
  */
  template <class T,
	    class Prop1, class Storage1, class Allocator1,
	    class Prop2, class Storage2, class Allocator2,
	    class Prop4, class Allocator4>
  void MltAdd(const T alpha,
	      const SeldonTranspose& TransA,
	      const Matrix<T, Prop1, Storage1, Allocator1>& A,
	      const SeldonTranspose& TransB,
	      const Matrix<T, Prop2, Storage2, Allocator2>& B,
	      const T beta,
	      Matrix<T, Prop4, PETScMPIAIJ, Allocator4>& C)
  {
    Matrix<T> tmp;
    Copy(C, tmp);
    MltAdd(alpha, TransA, A, TransB, B, beta, tmp);
    Copy(tmp, C);
  }


  // MLTADD //
  ////////////


  /////////
  // ADD //


  //! Adds two matrices.
  /*! It performs the operation \f$ B = \alpha A + B \f$ where \f$ \alpha \f$
    is a scalar, and \f$ A \f$ and \f$ B \f$ are matrices.
    \param[in] alpha scalar.
    \param[in] A matrix.
    \param[in,out] B matrix, result of the addition of \a B (on entry) and \a
    A times \a alpha.
  */
  template<class T, class Prop1, class Allocator1,
	   class Prop2, class Allocator2>
  void Add(const T& alpha,
	   const Matrix<T, Prop1, PETScMPIAIJ, Allocator1>& A,
	   Matrix<T, Prop2, PETScMPIAIJ, Allocator2>& B)
  {
    MatAXPY(B.GetPetscMatrix(), alpha,
            A.GetPetscMatrix(), DIFFERENT_NONZERO_PATTERN);
  }


  // ADD //
  /////////


  ///////////////
  // TRANSPOSE //


  //! Matrix transposition.
  template<class T, class Prop, class Allocator>
  void Transpose(Matrix<T, Prop, PETScMPIAIJ, Allocator>& A)
  {
    MatTranspose(A.GetPetscMatrix(), MAT_REUSE_MATRIX, &(A.GetPetscMatrix()));
  }


  //! Matrix transposition.
  template<class T, class Prop, class Allocator>
  void Transpose(Matrix<T, Prop, PETScMPIDense, Allocator>& A)
  {
    MatTranspose(A.GetPetscMatrix(), MAT_REUSE_MATRIX, &(A.GetPetscMatrix()));
  }


  // TRANSPOSE //
  ///////////////


  ///////////
  // GETLU //


  //! Returns the LU factorization of a matrix.
  /*! It factorizes the matrix \a A into \a L and \a U, so that \f$ A = L U
    \f$, \a L is a lower triangular matrix with ones on the diagonal, and \a U
    is an upper triangular matrix. On exit, the LU factorization is stored
    inside \a A: \a L in the lower part and \a U in the upper part. The
    diagonal elements are those of \a U. The diagonal elements of \a L are
    known to be ones.
    \param[in,out] A on entry, the matrix to be factorized; on exit, the LU
    factorization.
    \sa Seldon::SolveLU(const Matrix<T0, Prop0, Storage0, Allocator0>& M,
    Vector<T1, Storage1, Allocator1>& Y)
  */
  template <class T, class Prop0, class Allocator0>
  void GetLU(Matrix<T, Prop0, PETScMPIDense, Allocator0>& A)
  {
    int i, p, q, k;
    T temp, zero;
    SetComplexZero(zero);

    int ma = A.GetM();

#ifdef SELDON_CHECK_BOUNDS
    int na = A.GetN();
    if (na != ma)
      throw WrongDim("GetLU(A)", "The matrix must be squared.");
#endif


    for (i = 0; i < ma; i++)
      {
	for (p = i; p < ma; p++)
	  {
	    temp = zero;
	    for (k = 0; k < i; k++)
	      temp += A(p, k) * A(k, i);
	    A.SetBuffer(p, i , A(p, i) - temp);
	  }
        A.Flush();
	for (q = i+1; q < ma; q++)
	  {
	    temp = zero;
	    for (k = 0; k < i; k++)
	      temp += A(i, k) * A(k, q);
	    A.SetBuffer(i, q, (A(i,q ) - temp) / A(i, i));
	  }
        A.Flush();
      }
  }


  // GETLU //
  ///////////


} // namespace Seldon.


#endif
