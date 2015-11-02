// Copyright (C) 2001-2010 Vivien Mallet
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


class SparseMatrixTest: public testing::Test
{
protected:
  size_t m_;
  size_t n_;
  size_t Nelement_;
  int Nloop_;

public:
  void SetIdentity()
  {
    Matrix<double, General, ColSparse> A_col(m_, n_);
    Matrix<double, General, RowSparse> A_row(m_, n_);
    Matrix<double> A_full(m_, n_);

    A_col.SetIdentity();
    A_row.SetIdentity();
    A_full.SetIdentity();

    for (size_t i = 0; i < m_; i++)
      for (size_t j = 0; j < n_; j++)
        if (i == j)
        {
          ASSERT_TRUE(A_col(i, j) == double(1));
          ASSERT_TRUE(A_row(i, j) == double(1));
          ASSERT_TRUE(A_full(i, j) == double(1));
        }
        else
        {
          ASSERT_TRUE(A_col(i, j) == double(0));
          ASSERT_TRUE(A_row(i, j) == double(0));
          ASSERT_TRUE(A_full(i, j) == double(0));
        }

    // Testing the function applied to a non-empty matrix.
    A_row.FillRand(m_ + n_);
    A_full.FillRand();
    A_row.SetIdentity();
    A_full.SetIdentity();

    for (size_t i = 0; i < m_; i++)
      for (size_t j = 0; j < n_; j++)
        if (i == j)
        {
          ASSERT_TRUE(A_row(i, j) == double(1));
          ASSERT_TRUE(A_full(i, j) == double(1));
        }
        else
        {
          ASSERT_TRUE(A_row(i, j) == double(0));
          ASSERT_TRUE(A_full(i, j) == double(0));
        }
  }

  void GetRowCol()
  {
    size_t i, j;
    int k;
    for (k = 0; k < Nloop_; k++)
      {
        Matrix<double, General, RowSparse> A(m_, n_);
        A.FillRand(Nelement_);

        /*** To full vectors ***/
        Vector<double> row, column;
        // Rows.
        for (i = 0; i < m_; i++)
          {
            GetRow(A, i, row);
            for (j = 0; j < n_; j++)
              ASSERT_TRUE(A(i, j) == row(j));
          }
        // Columns.
        for (j = 0; j < n_; j++)
          {
            GetCol(A, j, column);
            for (i = 0; i < m_; i++)
              ASSERT_TRUE(A(i, j) == column(i));
          }

        /*** To sparse vectors ***/

        Vector<double, VectSparse> row_sparse, column_sparse;
        // Rows.
        for (i = 0; i < m_; i++)
          {
            GetRow(A, i, row_sparse);
            for (j = 0; j < n_; j++)
              ASSERT_TRUE(A(i, j) == row_sparse(j));
          }
        // Columns.
        for (j = 0; j < n_; j++)
          {
            GetCol(A, j, column_sparse);
            for (i = 0; i < m_; i++)
              ASSERT_TRUE(A(i, j) == column_sparse(i));
          }
      }
  }


  void Conversion()
  {
    srand(time(NULL));

    size_t i, j;
    double value;

    for (int k = 0; k < Nloop_; k++)
      {
        Matrix<double> A_full(m_, n_);
        A_full.Zero();

        Matrix<double, General, ArrayRowSparse> A_array(m_, n_);
        for (size_t l = 0; l < Nelement_; l++)
          {
            i = rand() % m_;
            j = rand() % n_;
            value = double(rand());
            A_array.AddInteraction(i, j, value);
            A_full(i, j) += value;
          }

        Matrix<double, General, RowSparse> A;

        Copy(A_array, A);

        for (size_t i = 0; i < m_; i++)
          for (size_t j = 0; j < n_; j++)
            {
              ASSERT_TRUE(A_full(i, j) == A(i, j));
              ASSERT_TRUE(A_full(i, j) == A_array(i, j));
            }

        Matrix<double, General, ColSparse> A_col;

        Copy(A, A_col);

        for (size_t i = 0; i < m_; i++)
          for (size_t j = 0; j < n_; j++)
            ASSERT_TRUE(A_full(i, j) == A_col(i, j));

        Vector<size_t> row_index, col_index;
        Vector<double> value;
        ConvertMatrix_to_Coordinates(A_array, row_index, col_index, value);

        ConvertMatrix_from_Coordinates(row_index, col_index, value, A);

        for (size_t i = 0; i < m_; i++)
          for (size_t j = 0; j < n_; j++)
            ASSERT_TRUE(A_full(i, j) == A(i, j));
      }
  }


  void Permutation()
  {
    size_t i, j;
    double value;

    for (int k = 0; k < Nloop_; k++)
      {
        Matrix<double> A_full(m_, n_);
        A_full.Zero();

        Matrix<double, General, ArrayRowSparse> A_array(m_, n_);
        for (size_t l = 0; l < Nelement_; l++)
          {
            i = rand() % m_;
            j = rand() % n_;
            value = double(rand());
            A_array.AddInteraction(i, j, value);
            A_full(i, j) += value;
          }

        Matrix<double, General, RowSparse> A;
        Copy(A_array, A);
        Matrix<double, General, ColSparse> A_col;
        Copy(A_array, A_col);

        Vector<size_t> row_permutation(m_);
        row_permutation.Fill();
        size_t tmp;
        for (size_t l = 0; l < m_; l++)
          {
            i = rand() % m_;
            j = rand() % m_;
            tmp = row_permutation(i);
            row_permutation(i) = row_permutation(j);
            row_permutation(j) = tmp;
          }

        Vector<size_t> col_permutation(n_);
        col_permutation.Fill();
        for (size_t l = 0; l < n_; l++)
          {
            i = rand() % n_;
            j = rand() % n_;
            tmp = col_permutation(i);
            col_permutation(i) = col_permutation(j);
            col_permutation(j) = tmp;
          }

        ApplyInversePermutation(A, row_permutation, col_permutation);
        ApplyInversePermutation(A_col, row_permutation, col_permutation);

        for (i = 0; i < m_; i++)
          for (j = 0; j < n_; j++)
            {
              ASSERT_TRUE(A_full(i, j) == A(row_permutation(i),
                                               col_permutation(j)));
              ASSERT_TRUE(A_full(i, j) == A_col(row_permutation(i),
                                                   col_permutation(j)));
            }
      }
  }


  void Transposition()
  {
    srand(time(NULL));

    size_t i, j;
    double value;

    for (int k = 0; k < Nloop_; k++)
      {
  {
    Matrix<double, General, ArrayRowSparse> A_array(m_, n_),
            A_array_t(m_, n_);
    for (size_t l = 0; l < Nelement_; l++)
      {
        i = rand() % m_;
        j = rand() % n_;
        value = double(rand());
        A_array.AddInteraction(i, j, value);
      }

    Copy(A_array, A_array_t);
    Transpose(A_array_t);

    Matrix<double, General, RowSparse> A(m_, n_);
    Copy(A_array, A);

    Transpose(A);

    for (size_t i = 0; i < m_; i++)
      for (size_t j = 0; j < n_; j++)
        ASSERT_TRUE(A_array_t(j, i) == A(j, i));

    Transpose(A);

    for (size_t i = 0; i < m_; i++)
      for (size_t j = 0; j < n_; j++)
        ASSERT_TRUE(A_array(i, j) == A(i, j));
  }

  {
    Matrix<double, General, ArrayRowSparse> A_array(m_, n_),
            A_array_t(m_, n_);
    for (size_t l = 0; l < Nelement_; l++)
      {
        i = rand() % m_;
        j = rand() % n_;
        value = double(rand());
        A_array.AddInteraction(i, j, value);
      }

    Copy(A_array, A_array_t);
    Transpose(A_array_t);

    Matrix<double, General, RowSparse> A(m_, n_);
    Copy(A_array, A);

    Transpose(A);

    for (size_t i = 0; i < m_; i++)
      for (size_t j = 0; j < n_; j++)
        ASSERT_TRUE(A_array_t(j, i) == A(j, i));

    Transpose(A);

    for (size_t i = 0; i < m_; i++)
      for (size_t j = 0; j < n_; j++)
        ASSERT_TRUE(A_array(i, j) == A(i, j));
  }

  {
    Matrix<double, General, ArrayRowSparse> A_array(m_, n_),
            A_array_t(m_, n_);
    for (size_t l = 0; l < Nelement_; l++)
      {
        i = rand() % m_;
        j = rand() % n_;
        value = double(rand());
        A_array.AddInteraction(i, j, value);
      }

    Copy(A_array, A_array_t);
    Transpose(A_array_t);

    Matrix<double, General, RowSparse, MallocObject<double> > A(m_, n_);
    Copy(A_array, A);

    Transpose(A);

    for (size_t i = 0; i < m_; i++)
      for (size_t j = 0; j < n_; j++)
        ASSERT_TRUE(A_array_t(j, i) == A(j, i));

    Transpose(A);

    for (size_t i = 0; i < m_; i++)
      for (size_t j = 0; j < n_; j++)
        ASSERT_TRUE(A_array(i, j) == A(i, j));
  }

  {
    Matrix<double, General, ArrayRowSparse> A_array(m_, n_),
            A_array_t(m_, n_);
    for (size_t l = 0; l < Nelement_; l++)
      {
        i = rand() % m_;
        j = rand() % n_;
        value = double(rand());
        A_array.AddInteraction(i, j, value);
      }

    Copy(A_array, A_array_t);
    Transpose(A_array_t);

    Matrix<double, General, RowSparse, NewAlloc<double> > A(m_, n_);
    Copy(A_array, A);

    Transpose(A);

    for (size_t i = 0; i < m_; i++)
      for (size_t j = 0; j < n_; j++)
        ASSERT_TRUE(A_array_t(j, i) == A(j, i));

    Transpose(A);

    for (size_t i = 0; i < m_; i++)
      for (size_t j = 0; j < n_; j++)
        ASSERT_TRUE(A_array(i, j) == A(i, j));
  }
      }
  }


  void SetRowcol()
  {
    srand(time(NULL));

    size_t i, j,  Nrow_sparse, Ncol_sparse;
    int k;
    for (k = 0; k < Nloop_; k++)
      {
  {
    Matrix<double, General, RowSparse> A(m_, n_);
    A.FillRand(Nelement_);
    for (i = 0; i < m_; i++)
      {
        Nrow_sparse =  rand() % n_ + 1;
        Vector<double, VectSparse> row_sparse(Nrow_sparse);
        for (size_t l = 0; l < Nrow_sparse; l++)
    row_sparse.Index(l) = rand() % n_;
        row_sparse.FillRand();
        row_sparse.Assemble();

        SetRow(row_sparse, i, A);

        for (j = 0; j < n_; j++)
    ASSERT_TRUE(A(i, j) == row_sparse(j));
      }
  }
  {
    Matrix<double, General, RowSparse> A(m_, n_);
    A.FillRand(Nelement_);
    for (j = 0; j < n_; j++)
      {
        Ncol_sparse =  rand() % m_ + 1;
        Vector<double, VectSparse> col_sparse(Ncol_sparse);
        for (size_t l = 0; l < Ncol_sparse; l++)
    col_sparse.Index(l) = rand() % m_;
        col_sparse.FillRand();
        col_sparse.Assemble();

        SetCol(col_sparse, j, A);

        for (i = 0; i < m_; i++)
    ASSERT_TRUE(A(i, j) == col_sparse(i));
      }
  }
      }
  }

};

  TEST_F(SparseMatrixTest, TestIdentity)
  {
    m_ = 25;
    n_ = 10;
    SetIdentity();

    m_ = 10;
    n_ = 25;
    SetIdentity();

    m_ = 20;
    n_ = 20;
    SetIdentity();
  }


  TEST_F(SparseMatrixTest, TestGetRowCol)
  {
    srand(time(NULL));

    m_ = 25;
    n_ = 10;
    Nelement_ = 30;
    Nloop_ = 10;
    GetRowCol();

    m_ = 10;
    n_ = 25;
    Nelement_ = 1;
    Nloop_ = 10;
    GetRowCol();

    m_ = 20;
    n_ = 20;
    Nelement_ = 300;
    Nloop_ = 10;
    GetRowCol();
  }


  TEST_F(SparseMatrixTest, TestConversion)
  {
    m_ = 25;
    n_ = 10;
    Nelement_ = 30;
    Nloop_ = 10;
    Conversion();

    m_ = 50;
    n_ = 60;
    Nelement_ = 5;
    Nloop_ = 10;
    Conversion();

    m_ = 50;
    n_ = 50;
    Nelement_ = 1;
    Nloop_ = 10;
    Conversion();

    m_ = 10;
    n_ = 25;
    Nelement_ = 30;
    Nloop_ = 10;
    Conversion();

    m_ = 10;
    n_ = 5;
    Nelement_ = 40;
    Nloop_ = 10;
    Conversion();

    m_ = 5;
    n_ = 10;
    Nelement_ = 100;
    Nloop_ = 2;
    Conversion();
  }


  TEST_F(SparseMatrixTest, TestPermutation)
  {
    m_ = 25;
    n_ = 10;
    Nelement_ = 30;
    Nloop_ = 10;
    Permutation();

    m_ = 50;
    n_ = 60;
    Nelement_ = 5;
    Nloop_ = 10;
    Permutation();

    m_ = 50;
    n_ = 50;
    Nelement_ = 1;
    Nloop_ = 10;
    Permutation();

    m_ = 10;
    n_ = 25;
    Nelement_ = 30;
    Nloop_ = 10;
    Permutation();

    m_ = 10;
    n_ = 5;
    Nelement_ = 40;
    Nloop_ = 10;
    Permutation();

    m_ = 5;
    n_ = 10;
    Nelement_ = 100;
    Nloop_ = 2;
    Permutation();
  }


  TEST_F(SparseMatrixTest, TestTransposition)
  {
    m_ = 25;
    n_ = 10;
    Nelement_ = 30;
    Nloop_ = 10;
    Transposition();

    m_ = 50;
    n_ = 60;
    Nelement_ = 5;
    Nloop_ = 10;
    Transposition();

    m_ = 50;
    n_ = 50;
    Nelement_ = 1;
    Nloop_ = 10;
    Transposition();

    m_ = 10;
    n_ = 25;
    Nelement_ = 30;
    Nloop_ = 10;
    Transposition();

    m_ = 10;
    n_ = 5;
    Nelement_ = 40;
    Nloop_ = 10;
    Transposition();

    m_ = 5;
    n_ = 10;
    Nelement_ = 100;
    Nloop_ = 2;
    Transposition();
  }


  TEST_F(SparseMatrixTest, TestSetRowCol)
  {
    m_ = 15;
    n_ = 5;
    Nelement_ = 20;
    Nloop_ = 10;
    SetRowcol();

    m_ = 50;
    n_ = 60;
    Nelement_ = 5;
    Nloop_ = 10;
    SetRowcol();

    m_ = 50;
    n_ = 50;
    Nelement_ = 1;
    Nloop_ = 10;
    SetRowcol();

    m_ = 10;
    n_ = 25;
    Nelement_ = 30;
    Nloop_ = 10;
    SetRowcol();

    m_ = 10;
    n_ = 5;
    Nelement_ = 40;
    Nloop_ = 10;
    SetRowcol();
  }


  


  

