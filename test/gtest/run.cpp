// Copyright (C) 2001-2010 Vivien Mallet, Marc Fragu
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

#define SELDON_WITH_BLAS
#define SELDON_WITH_LAPACK


#include "gtest/gtest.h"
#include "Seldon.hxx"
#include "matrix.hpp"
#include "array3d.hpp"
#include "lapack.hpp"
#include "sparse_matrix.hpp"
#include "sparse_linear_algebra.hpp"
#include "submatrix.hpp"
#include "vector_collection.hpp"
#include "heterogeneous_collection.hpp"


// Main function used to launch the gtest framework.
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
