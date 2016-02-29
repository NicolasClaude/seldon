// NewAlloc as default allocator in order to avoid problems
// with vectors of complex types
#define SELDON_DEFAULT_ALLOCATOR NewAlloc
// seldon will call abort() when encountering an exception
#define SELDON_WITH_ABORT
// no call of srand by Seldon
#define SELDON_WITHOUT_REINIT_RANDOM

#include "matrix_sparse/BandMatrix.hxx"
#include "matrix_sparse/BandMatrix.cxx"

using namespace Seldon;
typedef double Real_wp;
typedef complex<double> Complex_wp;

double threshold;

class MatrixBandTest: public testing::Test
{

public:
  static void SetUpTestCase()
  {
    threshold = 1e-4;
  }
};


template<class T>
void CheckSize(Matrix<T, General, BandedCol>& A, size_t n, size_t kl, size_t ku)
{
  A.Reallocate(n, n, kl, ku);
  ASSERT_EQ(n, A.GetM());
  ASSERT_EQ(n, A.GetN());
  ASSERT_EQ(kl, A.GetKL());
  ASSERT_EQ(ku, A.GetKU());
  ASSERT_EQ(n * (kl + ku + 1), A.GetDataSize());
}

template<class T>
void CheckClear(Matrix<T, General, BandedCol>& A, size_t n, size_t kl, size_t ku)
{
  A.Reallocate(n, n, kl, ku);
  A.Clear();
  ASSERT_EQ(0, A.GetM());
  ASSERT_EQ(0, A.GetDataSize());
}

template<class T>
void CheckZero(Matrix<T, General, BandedCol>& A, size_t n, size_t kl, size_t ku)
{
  A.Reallocate(n, n, kl, ku);
  A.Zero();
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      ASSERT_NEAR(abs(A(i, j)), 0, threshold);
}

template<class T>
void CheckGet(Matrix<T, General, BandedCol>& A, size_t n, size_t kl, size_t ku)
{
  A.Reallocate(n, n, kl, ku);
  typedef typename ClassComplexType<T>::Treal Treal;
  T zero; SetComplexZero(zero);


  Matrix<T, General, RowMajor> B(n, n);
  B.FillRand();
  Mlt(Treal(1e-9), B);

  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      {
	if ((i <= j+kl) && (i + ku >= j))
	  {
	    A.Get(i, j) = B(i, j);
	  }
	else
	  B(i, j) = zero;
      }

  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      ASSERT_EQ(A(i, j), B(i,j));
}

template<class T>
void CheckVal(Matrix<T, General, BandedCol>& A, size_t n, size_t kl, size_t ku)
{
  A.Reallocate(n, n, kl, ku);
  typedef typename ClassComplexType<T>::Treal Treal;
  T zero; SetComplexZero(zero);

  Matrix<T, General, RowMajor> B(n, n);
  B.FillRand();
  Mlt(Treal(1e-9), B);

  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      if ((i <= j+kl) && (i +ku >= j))
	A.Val(i, j) = B(i, j);
      else
        B(i, j) = 0;

  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      ASSERT_EQ(A(i, j), B(i,j));
}

template<class T>
void CheckAddInteraction(Matrix<T, General, BandedCol>& A, size_t n, size_t kl, size_t ku)
{
  A.Reallocate(n, n, kl, ku);
  typedef typename ClassComplexType<T>::Treal Treal;
  T zero; SetComplexZero(zero);

  Matrix<T, General, RowMajor> B(n, n);
  B.FillRand();
  Mlt(Treal(1e-9), B);

  Matrix<T, General, RowMajor> C(n, n);
  C.FillRand();
  Mlt(Treal(1e-9), C);

  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      if ((i <= j+kl) && (i + ku >= j))
	A.Val(i, j) = B(i, j);
      else
        B(i, j) = 0;


  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      {
	if ((i <= j+kl) && (i +ku  >= j))
	  A.AddInteraction(i, j, C(i, j));
	else
	  C(i, j) = zero;
      }

  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      ASSERT_EQ(A(i,j), B(i,j) + C(i,j));
}

template<class T>
void CheckAddInteractionRow(Matrix<T, General, BandedCol>& A, size_t n, size_t kl, size_t ku)
{
  A.Reallocate(n, n, kl, ku);
  typedef typename ClassComplexType<T>::Treal Treal;
  T zero; SetComplexZero(zero);

  Matrix<T, General, RowMajor> B(n, n);
  B.FillRand();
  Mlt(Treal(1e-9), B);

  Matrix<T, General, RowMajor> C(n, n);
  C.FillRand();
  Mlt(Treal(1e-9), C);

  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      if ((i <= j+kl) && (i + ku >= j))
	A.Val(i, j) = B(i, j);
      else
        B(i, j) = 0;

  IVect num(kl+ku+1); Vector<T> val(kl+ku+1);
  for (size_t i = 0; i < n; i++)
    {
      size_t nb_val = 0;
      for (size_t j = 0; j < n; j++)
        {
          if ((i <= j+kl) && (i +ku >= j ))
            {
              num(nb_val) = j;
              val(nb_val) = C(i, j);
              nb_val++;
            }
          else
            C(i, j) = 0;
        }
      A.AddInteractionRow(i, nb_val, num, val);
    }


  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      ASSERT_EQ(A(i,j), B(i,j) + C(i,j));
}

template<class T>
void CheckOperatorMlt(Matrix<T, General, BandedCol>& A, size_t n, size_t kl, size_t ku)
{
  A.Reallocate(n, n, kl, ku);
  typedef typename ClassComplexType<T>::Treal Treal;
  T zero; SetComplexZero(zero);
  Matrix<T, General, RowMajor> B(n, n);
  B.FillRand();
  Mlt(Treal(1e-9), B);

  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      if ((i <= j+kl) && (i + ku >= j))
         A.Val(i, j) = B(i, j);
      else
        B(i, j) = zero;

  T coef_mlt;
  SetComplexReal(to_num<Treal>("1.3"), coef_mlt);
  A *= coef_mlt;

  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
       ASSERT_EQ(A(i,j), B(i,j)* coef_mlt);
}

template<class T>
void CheckClearRow(Matrix<T, General, BandedCol>& A, size_t n, size_t kl, size_t ku)
{
  A.Reallocate(n, n, kl, ku);
  typedef typename ClassComplexType<T>::Treal Treal;
  T zero; SetComplexZero(zero);
  Matrix<T, General, RowMajor> B(n, n);
  B.FillRand();
  Mlt(Treal(1e-9), B);
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      if ((i <= j + kl) && (i + ku >= j))
         A.Val(i, j) = B(i, j);
      else
        B(i, j) = 0;

  A.ClearRow(5);
  for (size_t j = 0; j < n; j++)
    B(5, j) = zero;

  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
       ASSERT_EQ(A(i, j), B(i, j));

}

template<class T>
void CheckSetIdentity(Matrix<T, General, BandedCol>& A, size_t n, size_t kl, size_t ku)
{
  A.Reallocate(n, n, kl, ku);
  A.SetIdentity();
  T zero; SetComplexZero(zero);
  T one; SetComplexOne(one);
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      {
        T val_ex = zero;
        if (i == j)
          val_ex = one;
        ASSERT_EQ(A(i,j), val_ex);
      }
}

template<class T>
void CheckFill(Matrix<T, General, BandedCol>& A, size_t n, size_t kl, size_t ku)
{
  typedef typename ClassComplexType<T>::Treal Treal;
  A.Reallocate(n, n, kl, ku);
  T coef_alpha;
  SetComplexReal(to_num<Treal>("2.5"), coef_alpha);
  A.Fill(coef_alpha);
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      if ((i <= j+kl) && (i + ku >= j))
        ASSERT_EQ(A(i,j), coef_alpha);
}

TEST_F(MatrixBandTest, testGetSizes)
{

  Matrix<Real_wp, General, BandedCol> A;
  size_t n = 20, kl = 2, ku = 3;
  CheckSize(A, n, kl, ku);

  Matrix<Real_wp, General, BandedCol> B;
  n = 25;
  kl = 4;
  ku = 2;
  CheckSize(B, n, kl, ku);

  Matrix<Complex_wp, General, BandedCol> C;
  n = 20;
  kl = 2;
  ku = 3;
  CheckSize(C, n, kl, ku);

  Matrix<Complex_wp, General, BandedCol> D;
  n = 25;
  kl = 4;
  ku = 2;
  CheckSize(D, n, kl, ku);

}

TEST_F(MatrixBandTest, testClear)
{

  Matrix<Real_wp, General, BandedCol> A;
  size_t n = 20, kl = 2, ku = 3;
  CheckClear(A, n, kl, ku);

  Matrix<Real_wp, General, BandedCol> B;
  n = 25;
  kl = 4;
  ku = 2;
  CheckClear(B, n, kl, ku);

  Matrix<Complex_wp, General, BandedCol> C;
  n = 20;
  kl = 2;
  ku = 3;
  CheckClear(C, n, kl, ku);

  Matrix<Complex_wp, General, BandedCol> D;
  n = 25;
  kl = 4;
  ku = 2;
  CheckClear(D, n, kl, ku);

}

TEST_F(MatrixBandTest, testZero)
{

  Matrix<Real_wp, General, BandedCol> A;
  size_t n = 20, kl = 2, ku = 3;
  CheckZero(A, n, kl, ku);

  Matrix<Real_wp, General, BandedCol> B;
  n = 25;
  kl = 4;
  ku = 2;
  CheckZero(B, n, kl, ku);

  Matrix<Complex_wp, General, BandedCol> C;
  n = 20;
  kl = 2;
  ku = 3;
  CheckZero(C, n, kl, ku);

  Matrix<Complex_wp, General, BandedCol> D;
  n = 25;
  kl = 4;
  ku = 2;
  CheckZero(D, n, kl, ku);

}

TEST_F(MatrixBandTest, testGet)
{

  Matrix<Real_wp, General, BandedCol> A;
  size_t n = 20, kl = 2, ku = 3;
  CheckGet(A, n, kl, ku);

  Matrix<Real_wp, General, BandedCol> B;
  n = 25;
  kl = 4;
  ku = 2;
  CheckGet(B, n, kl, ku);

  Matrix<Complex_wp, General, BandedCol> C;
  n = 20;
  kl = 2;
  ku = 3;
  CheckGet(C, n, kl, ku);

  Matrix<Complex_wp, General, BandedCol> D;
  n = 25;
  kl = 4;
  ku = 2;
  CheckGet(D, n, kl, ku);

}

TEST_F(MatrixBandTest, testVal)
{

  Matrix<Real_wp, General, BandedCol> A;
  size_t n = 20, kl = 2, ku = 3;
  CheckVal(A, n, kl, ku);

  Matrix<Real_wp, General, BandedCol> B;
  n = 25;
  kl = 4;
  ku = 2;
  CheckVal(B, n, kl, ku);

  Matrix<Complex_wp, General, BandedCol> C;
  n = 20;
  kl = 2;
  ku = 3;
  CheckVal(C, n, kl, ku);

  Matrix<Complex_wp, General, BandedCol> D;
  n = 25;
  kl = 4;
  ku = 2;
  CheckVal(D, n, kl, ku);

}

TEST_F(MatrixBandTest, testAddInteraction)
{

  Matrix<Real_wp, General, BandedCol> A;
  size_t n = 20, kl = 2, ku = 3;
  CheckAddInteraction(A, n, kl, ku);

  Matrix<Real_wp, General, BandedCol> B;
  n = 25;
  kl = 4;
  ku = 2;
  CheckAddInteraction(B, n, kl, ku);

  Matrix<Complex_wp, General, BandedCol> C;
  n = 20;
  kl = 2;
  ku = 3;
  CheckAddInteraction(C, n, kl, ku);

  Matrix<Complex_wp, General, BandedCol> D;
  n = 25;
  kl = 4;
  ku = 2;
  CheckAddInteraction(D, n, kl, ku);

}

TEST_F(MatrixBandTest, testAddInteractionRow)
{

  Matrix<Real_wp, General, BandedCol> A;
  size_t n = 20, kl = 2, ku = 3;
  CheckAddInteractionRow(A, n, kl, ku);

  Matrix<Real_wp, General, BandedCol> B;
  n = 25;
  kl = 4;
  ku = 2;
  CheckAddInteractionRow(B, n, kl, ku);

  Matrix<Complex_wp, General, BandedCol> C;
  n = 20;
  kl = 2;
  ku = 3;
  CheckAddInteractionRow(C, n, kl, ku);

  Matrix<Complex_wp, General, BandedCol> D;
  n = 25;
  kl = 4;
  ku = 2;
  CheckAddInteractionRow(D, n, kl, ku);

}

TEST_F(MatrixBandTest, testOperatorMlt)
{

  Matrix<Real_wp, General, BandedCol> A;
  size_t n = 20, kl = 2, ku = 3;
  CheckOperatorMlt(A, n, kl, ku);

  Matrix<Real_wp, General, BandedCol> B;
  n = 25;
  kl = 4;
  ku = 2;
  CheckOperatorMlt(B, n, kl, ku);

  Matrix<Complex_wp, General, BandedCol> C;
  n = 20;
  kl = 2;
  ku = 3;
  CheckOperatorMlt(C, n, kl, ku);

  Matrix<Complex_wp, General, BandedCol> D;
  n = 25;
  kl = 4;
  ku = 2;
  CheckOperatorMlt(D, n, kl, ku);

}

TEST_F(MatrixBandTest, testClearRow)
{

  Matrix<Real_wp, General, BandedCol> A;
  size_t n = 20, kl = 2, ku = 3;
  CheckClearRow(A, n, kl, ku);

  Matrix<Real_wp, General, BandedCol> B;
  n = 25;
  kl = 4;
  ku = 2;
  CheckClearRow(B, n, kl, ku);

  Matrix<Complex_wp, General, BandedCol> C;
  n = 20;
  kl = 2;
  ku = 3;
  CheckClearRow(C, n, kl, ku);

  Matrix<Complex_wp, General, BandedCol> D;
  n = 25;
  kl = 4;
  ku = 2;
  CheckClearRow(D, n, kl, ku);

}


TEST_F(MatrixBandTest, testSetIdentity)
{

  Matrix<Real_wp, General, BandedCol> A;
  size_t n = 20, kl = 2, ku = 3;
  CheckSetIdentity(A, n, kl, ku);

  Matrix<Real_wp, General, BandedCol> B;
  n = 25;
  kl = 4;
  ku = 2;
  CheckSetIdentity(B, n, kl, ku);

  Matrix<Complex_wp, General, BandedCol> C;
  n = 20;
  kl = 2;
  ku = 3;
  CheckSetIdentity(C, n, kl, ku);

  Matrix<Complex_wp, General, BandedCol> D;
  n = 25;
  kl = 4;
  ku = 2;
  CheckSetIdentity(D, n, kl, ku);

}

TEST_F(MatrixBandTest, testFill)
{
  Matrix<Real_wp, General, BandedCol> A;
  size_t n = 20, kl = 2, ku = 3;
  CheckFill(A, n, kl, ku);

  Matrix<Real_wp, General, BandedCol> B;
  n = 25;
  kl = 4;
  ku = 2;
  CheckFill(B, n, kl, ku);

  Matrix<Complex_wp, General, BandedCol> C;
  n = 20;
  kl = 2;
  ku = 3;
  CheckFill(C, n, kl, ku);

  Matrix<Complex_wp, General, BandedCol> D;
  n = 25;
  kl = 4;
  ku = 2;
  CheckFill(D, n, kl, ku);

}
/*
TEST_F(MatrixBandTest, testPrint)
{
  Matrix<Real_wp, General, BandedCol> A;
  A.Reallocate(9,9,2,3);
  A.Fill(1);
  for (size_t i = 0; i < 9; i++)
    for (size_t j = 0; j < 9; j++)
      if ((i <= j+2) && (i + 3 >= j))
        A.Val(i, j) = i+j-9;
    for(unsigned int i = 0 ; i <9; i++)
    {
      for(size_t j=0; j < 9 ; j++)
        cout << i << " " << j << " " << A(i,j) << " | ";
      cout << endl ;
      }

}
*/
