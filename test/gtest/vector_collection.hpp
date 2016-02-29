// Copyright (C) 2001-2009 Marc Fragu
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


#include <cstdlib>
#include <ctime>

#include "Seldon.hxx"
#include "SeldonSolver.hxx"
#include "vector/VectorCollection.cxx"
using namespace Seldon;


class VectorCollectionTest: public testing::Test
{


protected:
  int Nloop_;
  size_t Nvector_;
  size_t Nsub_vector_max_;
  size_t m_;

public:
void mlt()
  {
    srand(time(NULL));

    size_t length;
    typedef double real;

    for (int N = 0; N < Nloop_; N++)
      {
	{  // For dense vector collection.
	  typedef Vector<real> vector_real_dense;

	  Vector<vector_real_dense, Collection> mlt;
	  Vector<vector_real_dense, Collection> A;
	  vector_real_dense U, W;
	  real alpha;
	  alpha = rand();
	  for (size_t k = 0; k < Nvector_; k++)
	    {
	      length = rand() % Nsub_vector_max_ + 1;

	      U.Reallocate(length);
	      U.FillRand();

	      A.AddVector(U);

	      W.Reallocate(length);
	      W.Copy(U);

	      Mlt(alpha, W);

	      mlt.AddVector(W);

	      U.Nullify();
	      W.Nullify();
	    }

	  Mlt(alpha, A);

	  for (size_t j = 0; j < Nvector_; j++)
	    for (size_t l = 0; l < A.GetVectorLength()(j); l++)
	      ASSERT_TRUE(A.GetVector(j)(l) == mlt.GetVector(j)(l));

	  for (size_t k = 0; k < Nvector_; k++)
	    A.GetVector()(k).~Vector();
	  for (size_t k = 0; k < Nvector_; k++)
	    mlt.GetVector()(k).~Vector();
	}

	{
	  // For sparse vector collection.

	  typedef Vector<real, VectSparse, NewAlloc<real> >
	    vector_real_sparse;

	  Vector<vector_real_sparse, Collection> mlt;
	  Vector<vector_real_sparse, Collection> A;
	  vector_real_sparse U, W;
	  real alpha;
	  alpha = rand();
	  for (size_t k = 0; k < Nvector_; k++)
	    {
	      length = rand() % Nsub_vector_max_ + 1;

	      U.Reallocate(length);
	      for (size_t l = 0; l < length; l++)
		{
		  U.Index(l) = rand() % length;
		  U.Value(l) = rand();
		}
	      U.Assemble();

	      W.Reallocate(length);
	      W.Copy(U);

	      A.AddVector(U);

	      Mlt(alpha, W);

	      mlt.AddVector(W);

	      U.Nullify();
	      W.Nullify();
	    }

	  Mlt(alpha, A);

	  for (size_t j = 0; j < Nvector_; j++)
	    for (size_t l = 0; l < A.GetVectorLength()(j); l++)
	      {
		ASSERT_TRUE(A.GetVector(j).Index(l) ==
			       mlt.GetVector(j).Index(l));
		ASSERT_TRUE(A.GetVector(j).Value(l) ==
			       mlt.GetVector(j).Value(l));
	      }

	  for (size_t k = 0; k < Nvector_; k++)
	    A.GetVector(k).~Vector();
	  for (size_t k = 0; k < Nvector_; k++)
	    mlt.GetVector(k).~Vector();
	}

      }
  }


  void add()
  {
    srand(time(NULL));

    size_t length;
    typedef double real;

    for (int N = 0; N < Nloop_; N++)
      {
	{  // Add two dense vector collections.
	  typedef Vector<real> vector_real_dense;

	  Vector<vector_real_dense, Collection> sum;
	  Vector<vector_real_dense, Collection> A, B;
	  vector_real_dense U, V, W, X;
	  for (size_t k = 0; k < Nvector_; k++)
	    {
	      length = rand() % Nsub_vector_max_ + 1;

	      U.Reallocate(length);
	      U.FillRand();
	      A.AddVector(U);

	      W.Reallocate(length);
	      W.Copy(U);

	      V.Reallocate(length);
	      V.FillRand();
	      B.AddVector(V);

	      X.Reallocate(length);
	      X.Copy(V);

	      Add(1.0, W, X);

	      sum.AddVector(X);

	      U.Nullify();
	      V.Nullify();
	      X.Nullify();
	    }

	  Add(1.0, A, B);

	  for (size_t j = 0; j < Nvector_; j++)
	    for (size_t l = 0; l < A.GetVectorLength()(j); l++)
	      ASSERT_TRUE(B.GetVector(j)(l) == sum.GetVector(j)(l));

	  for (size_t k = 0; k < Nvector_; k++)
	    A.GetVector()(k).~Vector();
	  for (size_t k = 0; k < Nvector_; k++)
	    B.GetVector()(k).~Vector();
	  for (size_t k = 0; k < Nvector_; k++)
	    sum.GetVector()(k).~Vector();
	}

	{
	  // Add two sparse vector collections.

	  typedef Vector<real, VectSparse, NewAlloc<real> >
	    vector_real_sparse;

	  Vector<vector_real_sparse, Collection> sum;
	  Vector<vector_real_sparse, Collection> A, B;
	  vector_real_sparse U, V, W, X;
	  for (size_t k = 0; k < Nvector_; k++)
	    {
	      length = rand() % Nsub_vector_max_ + 1;

	      U.Reallocate(length);
	      for (size_t l = 0; l < length; l++)
		{
		  U.Index(l) = rand() % length;
		  U.Value(l) = rand();
		}
	      U.Assemble();

	      W.Reallocate(length);
	      W.Copy(U);

	      V.Reallocate(length);
	      for (size_t l = 0; l < length; l++)
		{
		  V.Index(l) = rand() % length;
		  V.Value(l) = rand();
		}
	      V.Assemble();

	      X.Reallocate(length);
	      X.Copy(V);

	      A.AddVector(U);
	      B.AddVector(V);

	      Add(1.0, W, X);

	      sum.AddVector(X);

	      U.Nullify();
	      V.Nullify();
	      X.Nullify();
	    }

	  Add(1.0, A, B);

	  for (size_t j = 0; j < Nvector_; j++)
	    for (size_t l = 0; l < A.GetVectorLength()(j); l++)
	      {
		ASSERT_TRUE(B.GetVector(j).Index(l) ==
			       sum.GetVector(j).Index(l));
		ASSERT_TRUE(B.GetVector(j).Value(l) ==
			       sum.GetVector(j).Value(l));
	      }

	  for (size_t k = 0; k < Nvector_; k++)
	    A.GetVector(k).~Vector();
	  for (size_t k = 0; k < Nvector_; k++)
	    B.GetVector(k).~Vector();
	  for (size_t k = 0; k < Nvector_; k++)
            sum.GetVector(k).~Vector();
	}

	{  // Add a dense vector to a dense vector collections.
	  typedef Vector<real> vector_real_dense;

	  Vector<vector_real_dense, Collection> sum;
	  Vector<vector_real_dense, Collection> A;
	  vector_real_dense B, U, V, W, X;
	  for (size_t k = 0; k < Nvector_; k++)
	    {
	      length = rand() % Nsub_vector_max_ + 1;

	      U.Reallocate(length);
	      U.FillRand();
	      A.AddVector(U);

	      W.Reallocate(length);
	      W.Copy(U);

	      V.Reallocate(length);
	      V.FillRand();
	      B.PushBack(V);

	      X.Reallocate(length);
	      X.Copy(V);

	      Add(1.0, W, X);

	      sum.AddVector(X);

	      U.Nullify();
	      X.Nullify();
	    }

	  Add(1.0, B, A);

	  for (size_t j = 0; j < Nvector_; j++)
	    for (size_t l = 0; l < A.GetVectorLength()(j); l++)
	      ASSERT_TRUE(A.GetVector(j)(l) == sum.GetVector(j)(l));

	  for (size_t k = 0; k < Nvector_; k++)
	    A.GetVector()(k).~Vector();
	  for (size_t k = 0; k < Nvector_; k++)
	    sum.GetVector()(k).~Vector();
	}

      }
  }


  void copy()
  {
    srand(time(NULL));

    size_t length;
    typedef double real;

    for (int N = 0; N < Nloop_; N++)
      {
	{  // For dense vector collection.

	  typedef Vector<real> vector_real_dense;

	  Vector<vector_real_dense, Collection> A, B;
	  vector_real_dense U, V;
	  for (size_t k = 0; k < Nvector_; k++)
	    {
	      length = rand() % Nsub_vector_max_ + 1;

	      U.Reallocate(length);
	      U.FillRand();
	      A.AddVector(U);

	      V.Reallocate(length);
	      V.FillRand();
	      B.AddVector(V);

	      U.Nullify();
	      V.Nullify();
	    }

	  B.Deallocate();

	  Copy(A, B);

	  for (size_t j = 0; j < Nvector_; j++)
	    for (size_t l = 0; l < A.GetVectorLength()(j); l++)
	      ASSERT_TRUE(A.GetVector(j)(l) == B.GetVector(j)(l));

	  A.Deallocate();

	}

	{
	  // For sparse vector collection.

	  typedef Vector<real, VectSparse, NewAlloc<real> >
	    vector_real_sparse;

	  Vector<vector_real_sparse, Collection> A, B;
	  vector_real_sparse U, V;
	  for (size_t k = 0; k < Nvector_; k++)
	    {
	      length = rand() % Nsub_vector_max_ + 1;

	      U.Reallocate(length);
	      for (size_t l = 0; l < length; l++)
		{
		  U.Index(l) = rand() % length;
		  U.Value(l) = rand();
		}
	      U.Assemble();

	      V.Reallocate(length);
	      for (size_t l = 0; l < length; l++)
		{
		  V.Index(l) = rand() % length;
		  V.Value(l) = rand();
		}
	      V.Assemble();

	      A.AddVector(U);
	      B.AddVector(V);

	      U.Nullify();
	      V.Nullify();
	    }

	  B.Deallocate();

	  Copy(A, B);

	  for (size_t j = 0; j < Nvector_; j++)
	    for (size_t l = 0; l < A.GetVectorLength()(j); l++)
	      {
		ASSERT_TRUE(B.GetVector(j).Index(l) ==
			       A.GetVector(j).Index(l));
		ASSERT_TRUE(B.GetVector(j).Value(l) ==
			       A.GetVector(j).Value(l));
	      }

	  A.Deallocate();
	}
      }
  }


  void dot_product()
  {
    srand(time(NULL));

    size_t length;
    typedef double real;

    for (int N = 0; N < Nloop_; N++)
      {
	{  // For dense vector collection.

	  typedef Vector<real> vector_real_dense;

	  Vector<vector_real_dense, Collection> A, B;
	  vector_real_dense U, V;
	  real dot1(0.), dot2;
	  for (size_t k = 0; k < Nvector_; k++)
	    {
	      length = rand() % Nsub_vector_max_ + 1;

	      U.Reallocate(length);
	      U.FillRand();
	      A.AddVector(U);

	      V.Reallocate(length);
	      V.FillRand();
	      B.AddVector(V);

	      dot1 += DotProd(U, V);

	      U.Nullify();
	      V.Nullify();
	    }

	  dot2 = DotProd(A, B);

	  ASSERT_TRUE(dot1 == dot2);

	  for (size_t k = 0; k < Nvector_; k++)
	    A.GetVector()(k).~Vector();
	  for (size_t k = 0; k < Nvector_; k++)
	    B.GetVector()(k).~Vector();
	}

	{
	  // For sparse vector collection.

	  typedef Vector<real, VectSparse, NewAlloc<real> >
	    vector_real_sparse;

	  Vector<vector_real_sparse, Collection> A, B;
	  real dot1 = real(0.), dot2;
	  vector_real_sparse U, V;
	  for (size_t k = 0; k < Nvector_; k++)
	    {
	      length = rand() % Nsub_vector_max_ + 1;

	      U.Reallocate(length);
	      for (size_t l = 0; l < length; l++)
		{
		  U.Index(l) = rand() % length;
		  U.Value(l) = rand();
		}
	      U.Assemble();

	      V.Reallocate(length);
	      for (size_t l = 0; l < length; l++)
		{
		  V.Index(l) = rand() % length;
		  V.Value(l) = rand();
		}
	      V.Assemble();

	      A.AddVector(U);
	      B.AddVector(V);

	      dot1 += DotProd(U, V);

	      U.Nullify();
	      V.Nullify();
	    }

	  dot2 = DotProd(A, B);

	  ASSERT_TRUE(dot1 == dot2);

	  for (size_t k = 0; k < Nvector_; k++)
	    A.GetVector(k).~Vector();
	  for (size_t k = 0; k < Nvector_; k++)
	    B.GetVector(k).~Vector();
	}
      }
  }


  void mlt_add()
  {
    srand(time(NULL));

    size_t length;
    typedef double real;

    for (int N = 0; N < Nloop_; N++)
      {
	{  // Add two dense vector collections.
	  typedef Vector<size_t> vector_real_dense;

	  Matrix<real> M;
	  real alpha, beta;
	  alpha = rand();
	  beta = rand();

	  Vector<vector_real_dense, Collection> A;
	  vector_real_dense A_dense, U, W;
	  for (size_t k = 0; k < Nvector_; k++)
	    {
	      length = rand() % Nsub_vector_max_ + 1;

	      U.Reallocate(length);
	      U.FillRand();
	      A.AddVector(U);

	      W.Reallocate(length);
	      W.Copy(U);
	      A_dense.PushBack(U);

	      U.Nullify();
	    }
	  vector_real_dense Y1(m_), Y2(m_);
	  Y1.FillRand();
	  Copy(Y1, Y2);

	  M.Reallocate(m_, A.GetM());
	  M.FillRand();

	  MltAdd(alpha, M, A, beta, Y1);

	  MltAdd(alpha, M, A_dense, beta, Y2);

	  for (size_t l = 0; l < m_; l++)
	    ASSERT_TRUE(Y1(l) == Y2(l));

	  for (size_t k = 0; k < Nvector_; k++)
	    A.GetVector()(k).~Vector();
	}
      }
  }


  void write_read()
  {
    srand(time(NULL));

    size_t length;
    Vector<size_t, VectFull, MallocAlloc<size_t> > length_vector;
    typedef double real;

    for (int N = 0; N < Nloop_; N++)
      {
	typedef Vector<real> vector_real_dense;

	Vector<vector_real_dense, Collection> A, B;

	length_vector.Clear();
	vector_real_dense U;
	for (size_t k = 0; k < Nvector_; k++)
	  {
	    length = rand() % Nsub_vector_max_ + 1;
	    length_vector.PushBack(length);
	    U.Reallocate(length);
	    U.FillRand();
	    A.AddVector(U);
	    U.Nullify();
	  }

	A.Write("test.bin");

	B.Read("test.bin", length_vector);

	for (size_t j = 0; j < Nvector_; j++)
	  for (size_t l = 0; l < A.GetVectorLength()(j); l++)
	    ASSERT_TRUE(A.GetVector(j)(l) == B.GetVector(j)(l));

	A.Deallocate();
	B.Deallocate();
      }
  }


  void label()
  {
    typedef Vector<double> vector_real_dense;
    size_t length;

    srand(time(NULL));

    Vector<vector_real_dense, Collection> A;
    vector_real_dense vector1, vector2;

    length = rand() % Nsub_vector_max_ + 1;
    vector1.Reallocate(length);
    vector1.FillRand();
    A.AddVector(vector1, "vector1");
    length = rand() % Nsub_vector_max_ + 1;
    vector2.Reallocate(length);
    vector2.FillRand();
    A.AddVector(vector2, "vector2");
    vector_real_dense U, V;

    U.Copy(A.GetVector("vector1"));

    V.Copy(A.GetVector("vector2"));

    for (size_t l = 0; l < U.GetM(); l++)
      ASSERT_TRUE(vector1(l) == U(l));
    for (size_t l = 0; l < V.GetM(); l++)
      ASSERT_TRUE(vector2(l) == V(l));

    A.Nullify();
  }


  void collection()
  {
    typedef double real;
    typedef Vector<real> vector_real_dense;
    typedef Vector<vector_real_dense, Collection>
      collection_real_dense;

    size_t length;
    Vector<collection_real_dense, Collection> C;
    collection_real_dense A, B, D;
    vector_real_dense U;

    for (size_t k = 0; k < Nvector_; k++)
      {
	length = rand() % Nsub_vector_max_ + 1;
	U.Reallocate(length);
	U.FillRand();
	A.AddVector(U);
	U.Nullify();
      }

    for (size_t k = 0; k < Nvector_; k++)
      {
	length = rand() % Nsub_vector_max_ + 1;
	U.Reallocate(length);
	U.FillRand();
	B.AddVector(U);
	U.Nullify();
      }

    D.Copy(A,false);
    C.AddVector(A);
    C.AddVector(B);
    C.AddVector(D);
    // Checking adress.
    for (size_t i = 0; i < A.GetNvector(); i++)
      ASSERT_TRUE(&A.GetVector(i).GetData()[0] ==
		     &C.GetVector(0).GetVector(i).GetData()[0]);

    // Checking adress.
    for (size_t i = 0; i < B.GetNvector(); i++)
      ASSERT_TRUE(&B.GetVector(i).GetData()[0] ==
		     &C.GetVector(1).GetVector(i).GetData()[0]);

    A.Deallocate();
    B.Deallocate();
  }


  void CopyDuplicate()
  {
    size_t length;
    typedef double real;
    typedef Vector<real> vector_real_dense;
    Vector<vector_real_dense, Collection> A, B;

    vector_real_dense U;
    for (size_t k = 0; k < Nvector_; k++)
      {
        length = rand() % Nsub_vector_max_ + 1;
        U.Reallocate(length);
        U.FillRand();
        A.AddVector(U);
        U.Nullify();
      }

    B.Copy(A);
  for (size_t k = 0; k < Nvector_; k++)
    for (size_t i = 0; i < A.GetVector(k).GetM(); i++)
      A.GetVector(k)(i) = i;

  for (size_t k = 0; k < Nvector_; k++)
    for (size_t i = 0; i < A.GetVector(k).GetM(); i++)
      ASSERT_NE(A.GetVector(k)(i), B.GetVector(k)(i));
  }
};


  TEST_F(VectorCollectionTest, TestMlt)
  {
    Nloop_ = 10;
    Nvector_ = 1000;
    Nsub_vector_max_ = 10;
    mlt();

    Nloop_ = 10;
    Nvector_ = 10;
    Nsub_vector_max_ = 10000;
    mlt();

    Nloop_ = 10;
    Nvector_ = 100;
    Nsub_vector_max_ = 100;
    mlt();
  }


  TEST_F(VectorCollectionTest, TestAdd)
  {
    Nloop_ = 10;
    Nvector_ = 1000;
    Nsub_vector_max_ = 10;
    add();

    Nloop_ = 10;
    Nvector_ = 10;
    Nsub_vector_max_ = 10000;
    add();

    Nloop_ = 10;
    Nvector_ = 100;
    Nsub_vector_max_ = 100;
    add();
  }


  TEST_F(VectorCollectionTest, TestCopy)
  {
    Nloop_ = 1;
    Nvector_ = 10;
    Nsub_vector_max_ = 100;
    copy();

    Nloop_ = 10;
    Nvector_ = 10;
    Nsub_vector_max_ = 10000;
    copy();

    Nloop_ = 10;
    Nvector_ = 100;
    Nsub_vector_max_ = 100;
    copy();
  }


  TEST_F(VectorCollectionTest, TestDotProduct)
  {
    Nloop_ = 10;
    Nvector_ = 1000;
    Nsub_vector_max_ = 10;
    dot_product();

    Nloop_ = 10;
    Nvector_ = 10;
    Nsub_vector_max_ = 10000;
    dot_product();

    Nloop_ = 10;
    Nvector_ = 100;
    Nsub_vector_max_ = 100;
    dot_product();
  }


  TEST_F(VectorCollectionTest, TestMltAdd)
  {
    Nloop_ = 10;
    Nvector_ = 1000;
    Nsub_vector_max_ = 10;
    m_ = 50;
    mlt_add();

    Nloop_ = 10;
    Nvector_ = 10;
    Nsub_vector_max_ = 10000;
    m_ = 10;
    mlt_add();

    Nloop_ = 10;
    Nvector_ = 100;
    Nsub_vector_max_ = 100;
    m_ = 100;
    mlt_add();
  }


  TEST_F(VectorCollectionTest, TestWriteRead)
  {
    Nloop_ = 1;
    Nvector_ = 5;
    Nsub_vector_max_ = 10;
    m_ = 50;
    write_read();
  }


  TEST_F(VectorCollectionTest, TestLabel)
  {
    Nloop_ = 10;
    Nvector_ = 1000;
    Nsub_vector_max_ = 10;
    label();
  }


  TEST_F(VectorCollectionTest, TestCollection)
  {
    Nvector_ = 2;
    Nsub_vector_max_ = 4;
    collection();
  }

  TEST_F(VectorCollectionTest, TestCopyDuplicate)
  {
    Nvector_ = 2;
    Nsub_vector_max_ = 4;
    CopyDuplicate();
  }
