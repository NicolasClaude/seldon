// Copyright (C) 2001-2009 Vivien Mallet
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



#include "Seldon.hxx"
#include "SeldonSolver.hxx"
using namespace Seldon;

typedef complex<float> complexfloat;
typedef complex<double> complexdouble;


class Matrix_Test: public testing::Test
{
protected:
  size_t m_;
  size_t n_;

public:
  virtual void SetUp()
  {
    m_ = 25;
    n_ = 10;
  }
};


  TEST_F(Matrix_Test, TestConstructor)
  {
    {
      Matrix<@real_complex, General, @storage_rectangular> M@storage_rectangular_@real_complex(m_, n_);
      ASSERT_TRUE(M@storage_rectangular_@real_complex.GetM() == m_);
      ASSERT_TRUE(M@storage_rectangular_@real_complex.GetN() == n_);
    }
    {
      Matrix<@real_complex, General, @storage_rectangular_full> M@storage_rectangular_full_@real_complex(0, n_);
      ASSERT_TRUE(M@storage_rectangular_full_@real_complex.GetM() == 0);
      ASSERT_TRUE(M@storage_rectangular_full_@real_complex.GetN() == n_);
    }
    {
      Matrix<@real_complex, General, @storage_rectangular_full> M@storage_rectangular_full_@real_complex(m_, 0);
      ASSERT_TRUE(M@storage_rectangular_full_@real_complex.GetM() == m_);
      ASSERT_TRUE(M@storage_rectangular_full_@real_complex.GetN() == 0);
    }
    {
      Matrix<@real_complex, General, @storage_full_real_complex> M@storage_full_real_complex_@real_complex(0, 0);
      ASSERT_TRUE(M@storage_full_real_complex_@real_complex.GetM() == 0);
      ASSERT_TRUE(M@storage_full_real_complex_@real_complex.GetN() == 0);
      Matrix<@complex, General, @storage_full_complex> M@storage_full_complex_@complex(0, 0);
      ASSERT_TRUE(M@storage_full_complex_@complex.GetM() == 0);
      ASSERT_TRUE(M@storage_full_complex_@complex.GetN() == 0);
    }
    {
      Matrix<@real_complex, General, @storage_full_real_complex> M@storage_full_real_complex_@real_complex(m_, m_);
      ASSERT_TRUE(M@storage_full_real_complex_@real_complex.GetM() == m_);
      ASSERT_TRUE(M@storage_full_real_complex_@real_complex.GetN() == m_);
      Matrix<@complex, General, @storage_full_complex> M@storage_full_complex_@complex(m_, m_);
      ASSERT_TRUE(M@storage_full_complex_@complex.GetM() == m_);
      ASSERT_TRUE(M@storage_full_complex_@complex.GetN() == m_);
    }
  }


  TEST_F(Matrix_Test, TestReallocate)
  {
    Matrix<@real_complex, General, @storage_rectangular_full> M@storage_rectangular_full_@real_complex(m_, n_);

    M@storage_rectangular_full_@real_complex.Reallocate(0, 0);
    ASSERT_TRUE(M@storage_rectangular_full_@real_complex.GetM() == 0);
    ASSERT_TRUE(M@storage_rectangular_full_@real_complex.GetN() == 0);
    M@storage_rectangular_full_@real_complex.Reallocate(m_, n_);
    ASSERT_TRUE(M@storage_rectangular_full_@real_complex.GetM() == m_);
    ASSERT_TRUE(M@storage_rectangular_full_@real_complex.GetN() == n_);
    M@storage_rectangular_full_@real_complex.Reallocate(m_, n_);
    ASSERT_TRUE(M@storage_rectangular_full_@real_complex.GetM() == m_);
    ASSERT_TRUE(M@storage_rectangular_full_@real_complex.GetN() == n_);
    M@storage_rectangular_full_@real_complex.Reallocate(2 * m_, n_);
    ASSERT_TRUE(M@storage_rectangular_full_@real_complex.GetM() == 2 * m_);
    ASSERT_TRUE(M@storage_rectangular_full_@real_complex.GetN() == n_);
    M@storage_rectangular_full_@real_complex.Reallocate(2 * m_, 2 * n_);
    ASSERT_TRUE(M@storage_rectangular_full_@real_complex.GetM() == 2 * m_);
    ASSERT_TRUE(M@storage_rectangular_full_@real_complex.GetN() == 2 * n_);
    M@storage_rectangular_full_@real_complex.Reallocate(m_, 0);
    ASSERT_TRUE(M@storage_rectangular_full_@real_complex.GetM() == m_);
    ASSERT_TRUE(M@storage_rectangular_full_@real_complex.GetN() == 0);
    M@storage_rectangular_full_@real_complex.Reallocate(0, n_);
    ASSERT_TRUE(M@storage_rectangular_full_@real_complex.GetM() == 0);
    ASSERT_TRUE(M@storage_rectangular_full_@real_complex.GetN() == n_);
    M@storage_rectangular_full_@real_complex.Reallocate(0, 0);
    ASSERT_TRUE(M@storage_rectangular_full_@real_complex.GetM() == 0);
    ASSERT_TRUE(M@storage_rectangular_full_@real_complex.GetN() == 0);
    M@storage_rectangular_full_@real_complex.Reallocate(0, n_);
    ASSERT_TRUE(M@storage_rectangular_full_@real_complex.GetM() == 0);
    ASSERT_TRUE(M@storage_rectangular_full_@real_complex.GetN() == n_);
    M@storage_rectangular_full_@real_complex.Reallocate(m_, n_);
    ASSERT_TRUE(M@storage_rectangular_full_@real_complex.GetM() == m_);
    ASSERT_TRUE(M@storage_rectangular_full_@real_complex.GetN() == n_);
    M@storage_rectangular_full_@real_complex.Reallocate(0, n_);
    ASSERT_TRUE(M@storage_rectangular_full_@real_complex.GetM() == 0);
    ASSERT_TRUE(M@storage_rectangular_full_@real_complex.GetN() == n_);
  }
