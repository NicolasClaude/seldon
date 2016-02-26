#ifndef SELDON_FILE_FUNCTIONS_MATVECT_PETSC_CXX
#define SELDON_FILE_FUNCTIONS_MATVECT_PETSC_CXX


namespace Seldon
{


  ////////////
  /// MLT  ///


  template <class T1, class Prop1, class Allocator1,
	    class T2, class Allocator2,
	    class T4, class Allocator4>
  void MltVector(const Matrix<T1, Prop1, PETScMPIDense, Allocator1>& M,
		 const Vector<T2, PETScPar, Allocator2>& X,
		 Vector<T4, PETScPar, Allocator4>& Y)
  {
    MatMult(M.GetPetscMatrix(), X.GetPetscVector(), Y.GetPetscVector());
  }


  /// MLT  ///
  ////////////


  ////////////
  // MLTADD //


  template <class T0,
            class T1, class Prop1, class Allocator1,
            class T2, class Allocator2,
            class T3,
            class T4, class Allocator4>
  void MltAddVector(const T0 alpha,
              const Matrix<T1, Prop1, PETScMPIAIJ, Allocator1>& M,
              const Vector<T2, PETScSeq, Allocator2>& X,
              const T3 beta, Vector<T4, PETScSeq, Allocator4>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(M, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif
    if (beta == T3(0))
      if (alpha == T0(0))
        {
          Y.Fill(T4(0));
          return;
        }
      else
        {
          MatMult(M.GetPetscMatrix(), X.GetPetscVector(), Y.GetPetscVector());
          if (alpha != T0(1))
            VecScale(Y.GetPetscVector(), alpha);
          return;
        }
    if (alpha == T0(1))
      {
        if (beta != T3(1))
          VecScale(Y.GetPetscVector(), beta);
        MatMultAdd(M.GetPetscMatrix(), X.GetPetscVector(),
                   Y.GetPetscVector(), Y.GetPetscVector());
        return;
      }
    Vector<T2, PETScSeq, Allocator2> tmp;
    MatMult(M.GetPetscMatrix(), X.GetPetscVector(), tmp.GetPetscVector());
    VecAXPBY(Y.GetPetscVector(), alpha, beta, tmp.GetPetscVector());
    return;
  }

  /*! \brief Performs the multiplication of a matrix with a vector, and adds
    the result to another vector.
   */
  /*! It performs the operation \f$ Y = \alpha M X + \beta Y \f$ where \f$
    \alpha \f$ and \f$ \beta \f$ are scalars, \f$ M \f$ is a \f$ m \times n
    \f$ matrix, and \f$ X \f$ is a vector of length \f$ n \f$. The vector \f$
    Y \f$ must be of length \f$ m \f$.
    \param[in] alpha scalar.
    \param[in] M m by n matrix.
    \param[in] X vector of length n.
    \param[in] beta scalar.
    \param[in,out] Y vector of length m, result of the product of \a M by \a
    X, times \a alpha, plus \a Y (on entry) times \a beta.
  */
  template <class T0,
	    class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Allocator4>
  void MltAddVector(const T0 alpha,
	      const Matrix<T1, Prop1, Storage1, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta,
	      Vector<T4, PETScPar, Allocator4>& Y)
  {
    int ma = M.GetM();
    int na = M.GetN();

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(M, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    Mlt(beta, Y);
    int l,h;
    VecGetOwnershipRange(Y.GetPetscVector(), &l, &h);
    T4 zero(0);
    T4 temp;
    T4 alpha_(alpha);
    for (int i = l; i < h; i++)
      {
	temp = zero;
	for (int j = 0; j < na; j++)
	  temp += M(i, j) * X(j);
        Y.SetBuffer(i, temp * alpha_ + Y(i));
      }
    Y.Flush();
  }


  template <class T0,
            class T1, class Prop1, class Allocator1,
            class T2, class Allocator2,
            class T3,
            class T4, class Allocator4>
  void MltAddVector(const T0 alpha,
              const Matrix<T1, Prop1, PETScMPIDense, Allocator1>& M,
              const Vector<T2, PETScPar, Allocator2>& X,
              const T3 beta, Vector<T4, PETScPar, Allocator4>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(M, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif
    if (beta == T3(0))
      if (alpha == T0(0))
        {
          Y.Fill(T4(0));
          return;
        }
      else
        {
          MatMult(M.GetPetscMatrix(), X.GetPetscVector(), Y.GetPetscVector());
          if (alpha != T0(1))
            VecScale(Y.GetPetscVector(), alpha);
          return;
        }
    if (alpha == T0(1))
      {
        if (beta != T3(1))
          VecScale(Y.GetPetscVector(), beta);
        MatMultAdd(M.GetPetscMatrix(), X.GetPetscVector(),
                   Y.GetPetscVector(),Y.GetPetscVector());
        return;
      }
    Vector<T2, PETScPar, Allocator2> tmp(X.GetM());

    MatMult(M.GetPetscMatrix(), X.GetPetscVector(), tmp.GetPetscVector());
    VecAXPBY(Y.GetPetscVector(), alpha, beta, tmp.GetPetscVector());
    return;
  }


  template <class T,
            class Prop1, class Allocator1,
            class Allocator2,
            class Allocator4>
  void MltAddVector(const T alpha,
              const Matrix<T, Prop1, PETScMPIAIJ, Allocator1>& M,
              const Vector<T, PETScPar, Allocator2>& X,
              const T beta, Vector<T, PETScPar, Allocator4>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(M, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif
    if (beta == T(0))
      if (alpha == T(0))
        {
          Y.Fill(T(0));
          return;
        }
      else
        {
          MatMult(M.GetPetscMatrix(), X.GetPetscVector(), Y.GetPetscVector());
          if (alpha != T(1))
            VecScale(Y.GetPetscVector(), alpha);
          return;
        }
    if (alpha == T(1))
      {
        if (beta != T(1))
          VecScale(Y.GetPetscVector(), beta);
        MatMultAdd(M.GetPetscMatrix(), X.GetPetscVector(),
                   Y.GetPetscVector(),Y.GetPetscVector());
        return;
      }
    Vector<T, PETScPar, Allocator2> tmp;
    tmp.Copy(Y);
    MatMult(M.GetPetscMatrix(), X.GetPetscVector(), tmp.GetPetscVector());
    VecAXPBY(Y.GetPetscVector(), alpha, beta, tmp.GetPetscVector());
    return;
  }


  template <class T0,
            class T1, class Prop1, class Allocator1,
            class T2, class Allocator2,
            class T3,
            class T4, class Allocator4>
  void MltAddVector(const T0 alpha,
              const Matrix<T1, Prop1, PETScMPIAIJ, Allocator1>& M,
              const Vector<T2, VectFull, Allocator2>& X,
              const T3 beta, Vector<T4, PETScSeq, Allocator4>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(M, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    Vector<T4, PETScSeq, Allocator4> X_Petsc;
    X_Petsc.Reallocate(X.GetM());
    for (int i = 0; i < X.GetM(); i++)
      X_Petsc.SetBuffer(i, X(i));
    X_Petsc.Flush();

    if (beta == T3(0))
      if (alpha == T0(0))
        {
          Y.Fill(T4(0));
          return;
        }
      else
        {
          MatMult(M.GetPetscMatrix(), X_Petsc.GetPetscVector(),
                  Y.GetPetscVector());
          if (alpha != T0(1))
            VecScale(Y.GetPetscVector(), alpha);
          return;
        }
    if (alpha == T0(1))
      {
        if (beta != T3(1))
          VecScale(Y.GetPetscVector(), beta);
        MatMultAdd(M.GetPetscMatrix(), X_Petsc.GetPetscVector(),
                   Y.GetPetscVector(),Y.GetPetscVector());
        return;
      }
    Vector<T2, PETScSeq, Allocator2> tmp;
    tmp.Copy(Y);
    MatMult(M.GetPetscMatrix(), X_Petsc.GetPetscVector(),
            tmp.GetPetscVector());
    VecAXPBY(Y.GetPetscVector(), alpha, beta, tmp.GetPetscVector());
    return;
  }


  template <class T0,
            class T1, class Prop1, class Allocator1,
            class T2, class Allocator2,
            class T3,
            class T4, class Allocator4>
  void MltAddVector(const T0 alpha,
              const Matrix<T0, Prop1, PETScMPIAIJ, Allocator1>& M,
              const Vector<T0, VectFull, Allocator2>& X,
              const T0 beta, Vector<T4, PETScPar, Allocator4>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(M, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    Vector<T4, PETScPar, Allocator4> X_Petsc;
    X_Petsc.SetCommunicator(M.GetCommunicator());
    X_Petsc.Reallocate(X.GetM());
    for (int i = 0; i < X.GetM(); i++)
      X_Petsc.SetBuffer(i, X(i));
    X_Petsc.Flush();

    if (beta == T3(0))
      if (alpha == T0(0))
        {
          Y.Fill(T4(0));
          return;
        }
      else
        {
          MatMult(M.GetPetscMatrix(), X_Petsc.GetPetscVector(),
                  Y.GetPetscVector());
          if (alpha != T0(1))
            VecScale(Y.GetPetscVector(), alpha);
          return;
        }
    if (alpha == T0(1))
      {
        if (beta != T3(1))
          VecScale(Y.GetPetscVector(), beta);
        MatMultAdd(M.GetPetscMatrix(), X_Petsc.GetPetscVector(),
                   Y.GetPetscVector(),Y.GetPetscVector());
        return;
      }
    Vector<T2, PETScPar, Allocator2> tmp;
    tmp.Copy(Y);
    MatMult(M.GetPetscMatrix(), X_Petsc.GetPetscVector(),
            tmp.GetPetscVector());
    VecAXPBY(Y.GetPetscVector(), alpha, beta, tmp.GetPetscVector());
    return;
  }


  template <class T0,
            class T1, class Prop1, class Allocator1,
            class T2, class Allocator2,
            class T3,
            class T4, class Allocator4>
  void MltAddVector(const T0 alpha,
              const Matrix<T1, Prop1, PETScMPIDense, Allocator1>& M,
              const Vector<T2, VectFull, Allocator2>& X,
              const T3 beta, Vector<T4, PETScPar, Allocator4>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(M, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    Vector<T4, PETScPar, Allocator4> X_Petsc;
     X_Petsc.SetCommunicator(M.GetCommunicator());
    X_Petsc.Reallocate(X.GetM());
    for (int i = 0; i < X.GetM(); i++)
      X_Petsc.SetBuffer(i, X(i));
    X_Petsc.Flush();

    if (beta == T3(0))
      if (alpha == T0(0))
        {
          Y.Fill(T4(0));
          return;
        }
      else
        {
          MatMult(M.GetPetscMatrix(), X_Petsc.GetPetscVector(),
                  Y.GetPetscVector());
          if (alpha != T0(1))
            VecScale(Y.GetPetscVector(), alpha);
          return;
        }
    if (alpha == T0(1))
      {
        if (beta != T3(1))
          VecScale(Y.GetPetscVector(), beta);
        MatMultAdd(M.GetPetscMatrix(), X_Petsc.GetPetscVector(),
                   Y.GetPetscVector(),Y.GetPetscVector());
        return;
      }
    Vector<T2, PETScPar, Allocator2> tmp;
    tmp.Copy(Y);
    tmp.SetCommunicator(M.GetCommunicator());
    MatMult(M.GetPetscMatrix(), X_Petsc.GetPetscVector(),
            tmp.GetPetscVector());
    VecAXPBY(Y.GetPetscVector(), alpha, beta, tmp.GetPetscVector());
    return;
  }


  template <class T0,
            class T1, class Prop1, class Allocator1,
            class T2, class Allocator2,
            class T3,
            class T4, class Allocator4>
  void MltAddVector(const T0 alpha,
              const Matrix<T1, Prop1, PETScMPIAIJ, Allocator1>& M,
              const Vector<T2, PETScPar, Allocator2>& X,
              const T3 beta, Vector<T4, VectFull, Allocator4>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(M, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    Vector<T4, PETScPar, Allocator4> Y_Petsc;
     Y_Petsc.SetCommunicator(M.GetCommunicator());
    Y_Petsc.Reallocate(Y.GetM());
    for (int i = 0; i < Y.GetM(); i++)
      Y_Petsc.SetBuffer(i, Y(i));
    Y_Petsc.Flush();
    MltAddVector(alpha, M, X, beta, Y_Petsc);
    for (int i = 0; i < Y.GetM(); i++)
      Y(i) = Y_Petsc(i);
    return;
  }


  /*! \brief Performs the multiplication of a matrix (possibly transposed)
    with a vector, and adds the result to another vector.
   */
  /*! It performs the operation \f$ Y = \alpha M X + \beta Y \f$ or \f$ Y =
    \alpha M^T X + \beta Y \f$ where \f$ \alpha \f$ and \f$ \beta \f$ are
    scalars, \f$ M \f$ is a \f$ m \times n \f$ matrix or a \f$ n \times m \f$
    matrix, and \f$ X \f$ is a vector of length \f$ n \f$. The vector \f$ Y
    \f$ must be of length \f$ m \f$.
    \param[in] alpha scalar.
    \param[in] Trans transposition status of \a M: it may be SeldonNoTrans or
    SeldonTrans.
    \param[in] M m by n matrix, or n by m matrix if transposed.
    \param[in] X vector of length n.
    \param[in] beta scalar.
    \param[in,out] Y vector of length m, result of the product of \a M by \a
    X, times \a alpha, plus \a Y (on entry) times \a beta.
    \note \a Trans must be either SeldonNoTrans or SeldonTrans: ConjTrans is
    not supported.
  */
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Allocator2,
	    class T3,
	    class T4, class Allocator4>
  void MltAddVector(const T0 alpha,
	      const SeldonTranspose& Trans,
	      const Matrix<T1, Prop1, PETScMPIAIJ, Allocator1>& M,
	      const Vector<T2, PETScPar, Allocator2>& X,
	      const T3 beta,
	      Vector<T4, PETScPar, Allocator4>& Y)
  {
    if (Trans.NoTrans())
      {
        MltAdd(alpha, M, X, beta, Y);
        return;
      }
    else if (Trans.ConjTrans())
      throw WrongArgument("MltAdd(alpha, trans, M, X, beta, Y)",
                          "Complex conjugation not supported.");

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Trans, M, X, Y, "MltAdd(alpha, trans, M, X, beta, Y)");
#endif
    Matrix<T1, Prop1, PETScMPIAIJ, Allocator1> temp;
    temp.Copy(M);
    Transpose(temp);
    MltAdd(alpha, temp, X, beta, Y);
  }


  /*! \brief Performs the multiplication of a matrix (possibly transposed)
    with a vector, and adds the result to another vector.
   */
  /*! It performs the operation \f$ Y = \alpha M X + \beta Y \f$ or \f$ Y =
    \alpha M^T X + \beta Y \f$ where \f$ \alpha \f$ and \f$ \beta \f$ are
    scalars, \f$ M \f$ is a \f$ m \times n \f$ matrix or a \f$ n \times m \f$
    matrix, and \f$ X \f$ is a vector of length \f$ n \f$. The vector \f$ Y
    \f$ must be of length \f$ m \f$.
    \param[in] alpha scalar.
    \param[in] Trans transposition status of \a M: it may be SeldonNoTrans or
    SeldonTrans.
    \param[in] M m by n matrix, or n by m matrix if transposed.
    \param[in] X vector of length n.
    \param[in] beta scalar.
    \param[in,out] Y vector of length m, result of the product of \a M by \a
    X, times \a alpha, plus \a Y (on entry) times \a beta.
    \note \a Trans must be either SeldonNoTrans or SeldonTrans: ConjTrans is
    not supported.
  */
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Allocator2,
	    class T3,
	    class T4, class Allocator4>
  void MltAddVector(const T0 alpha,
	      const SeldonTranspose& Trans,
	      const Matrix<T1, Prop1, PETScMPIDense, Allocator1>& M,
	      const Vector<T2, PETScPar, Allocator2>& X,
	      const T3 beta,
	      Vector<T4, PETScPar, Allocator4>& Y)
  {
    if (Trans.NoTrans())
      {
        MltAdd(alpha, M, X, beta, Y);
        return;
      }
    else if (Trans.ConjTrans())
      throw WrongArgument("MltAdd(alpha, trans, M, X, beta, Y)",
                          "Complex conjugation not supported.");

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Trans, M, X, Y, "MltAdd(alpha, trans, M, X, beta, Y)");
#endif
    Matrix<T1, Prop1, PETScMPIDense, Allocator1> tmp;
    tmp.Copy(M);
    Transpose(tmp);
    MltAdd(alpha, tmp, X, beta, Y);
  }


  /*! \brief Performs the multiplication of a matrix (possibly transposed)
    with a vector, and adds the result to another vector.
   */
  /*! It performs the operation \f$ Y = \alpha M X + \beta Y \f$ or \f$ Y =
    \alpha M^T X + \beta Y \f$ where \f$ \alpha \f$ and \f$ \beta \f$ are
    scalars, \f$ M \f$ is a \f$ m \times n \f$ matrix or a \f$ n \times m \f$
    matrix, and \f$ X \f$ is a vector of length \f$ n \f$. The vector \f$ Y
    \f$ must be of length \f$ m \f$.
    \param[in] alpha scalar.
    \param[in] Trans transposition status of \a M: it may be SeldonNoTrans or
    SeldonTrans.
    \param[in] M m by n matrix, or n by m matrix if transposed.
    \param[in] X vector of length n.
    \param[in] beta scalar.
    \param[in,out] Y vector of length m, result of the product of \a M by \a
    X, times \a alpha, plus \a Y (on entry) times \a beta.
    \note \a Trans must be either SeldonNoTrans or SeldonTrans: ConjTrans is
    not supported.
  */
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Allocator2,
	    class T3,
	    class T4, class Allocator4>
  void MltAddVector(const T0 alpha,
	      const SeldonTranspose& Trans,
	      const Matrix<T1, Prop1, PETScMPIAIJ, Allocator1>& M,
	      const Vector<T2, PETScPar, Allocator2>& X,
	      const T3 beta,
	      Vector<T4, VectFull, Allocator4>& Y)
  {
    if (Trans.NoTrans())
      {
        MltAdd(alpha, M, X, beta, Y);
        return;
      }
    else if (Trans.ConjTrans())
      throw WrongArgument("MltAdd(alpha, trans, M, X, beta, Y)",
                          "Complex conjugation not supported.");

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Trans, M, X, Y, "MltAdd(alpha, trans, M, X, beta, Y)");
#endif
    Matrix<T1, Prop1, PETScMPIAIJ, Allocator1> tmp;
    tmp.Copy(M);
    Transpose(tmp);
    MltAdd(alpha, tmp, X, beta, Y);
  }


  // MLTADD //
  ////////////


  /////////////
  // SOLVELU //


  //! Solves a linear system using LU factorization.
  /*! This function solves \f$ M X = Y \f$ where \f$ M \f$ is a matrix, and
    \f$ X \f$ and \f$ Y \f$ are vectors.
    \param[in] M the matrix of the linear system, to be factorized in LU
    form. On exit, \a M contains its LU factorization.
    \param[in,out] Y on entry, the right-hand side \f$ Y \f$; on exit, the
    solution \f$ X \f$ of the system.
  */
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Allocator1>
  void GetAndSolveLU(Matrix<T0, Prop0, PETScMPIDense, Allocator0>& M,
                     Vector<T1, PETScPar, Allocator1>& Y)
  {
    GetLU(M);
    SolveLU(M, Y);
  }


  //! Solves a linear system whose matrix has been LU-factorized.
  /*! This function solves \f$ M X = Y \f$ where \f$ M \f$ is a matrix, and
    \f$ X \f$ and \f$ Y \f$ are vectors. The matrix \a M cannot be provided as
    such to this function: it must already be factorized in LU form.
    \param[in] M the matrix of the linear system, already factorized in LU
    form. The lower part of \a M should be \a L, and the upper part should be
    \a U. The diagonal of \a M should be the diagonal of \a U. The diagonal
    elements of \a L are assumed to be ones.
    \param[in,out] Y on entry, the right-hand side \f$ Y \f$; on exit, the
    solution \f$ X \f$ of the system.
    \sa Seldon::GetLU(Matrix<T0, Prop0, Storage0, Allocator0>& A) to factorize
    a matrix before using this function.
  */
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Allocator1>
  void SolveLuVector(const Matrix<T0, Prop0, PETScMPIDense, Allocator0>& M,
		     Vector<T1, PETScPar, Allocator1>& Y)
  {
    int i, k;
    T1 tmp;

    int ma = M.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
    int na = M.GetN();
    if (na != ma)
      throw WrongDim("SolveLU(M, Y)",
		     "The matrix must be squared.");

    CheckDim(M, Y, "SolveLU(M, Y)");
#endif


    // Forward substitution.
    for (i = 0; i < ma; i++)
      {
	SetComplexZero(tmp);
	for (k = 0; k < i; k++)
	  tmp += M(i, k) * Y(k);
	Y.SetBuffer(i, (Y(i) - tmp) / M(i, i));
      }
    Y.Flush();

    // Back substitution.
    for (i = ma - 2; i > -1; i--)
      {
	SetComplexZero(tmp);
	for (k = i + 1; k < ma; k++)
	  tmp += M(i, k) * Y(k);
	Y.SetBuffer(i, Y(i) - tmp);
      }
    Y.Flush();
  }


  // SOLVELU //
  /////////////


}

#endif
