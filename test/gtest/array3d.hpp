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




#include "../../Seldon.hxx"
using namespace Seldon;


class Array3DTest: public testing::Test
  {



protected:
  size_t length1_;
  size_t length2_;
  size_t length3_;

public:
  void SetUp()
  {
    length1_ = 10;
    length2_ = 25;
    length3_ = 9;
  }

};

  TEST_F(Array3DTest, TestFill)
  {
    size_t i, j, k;
    Array3D<double> A(length1_, length2_, length3_);
    A.Fill();
    for (i = 0; i < length1_; i++)
      for (j = 0; j < length2_; j++)
	for (k = 0; k < length3_; k++)
	  ASSERT_TRUE(A(i, j, k) == double(i * length2_ * length3_
					      + j * length3_ + k));
  }


  TEST_F(Array3DTest, TestReallocate)
  {
    Array3D<double> A(length1_, length2_, length3_);
    A.Reallocate(length2_, length1_, length3_);
    ASSERT_TRUE(A.GetLength1() == length2_);
    ASSERT_TRUE(A.GetLength2() == length1_);
    ASSERT_TRUE(A.GetLength3() == length3_);
    A.Reallocate(0, length2_, 0);
    ASSERT_TRUE(A.GetLength1() == 0);
    ASSERT_TRUE(A.GetLength2() == length2_);
    ASSERT_TRUE(A.GetLength3() == 0);
    A.Reallocate(0, 0, 0);
    ASSERT_TRUE(A.GetLength1() == 0);
    ASSERT_TRUE(A.GetLength2() == 0);
    ASSERT_TRUE(A.GetLength3() == 0);
  }

