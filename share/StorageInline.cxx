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


#ifndef SELDON_FILE_STORAGE_INLINE_CXX

#include "Storage.hxx"

namespace Seldon
{


  //////////////////////
  // GENERAL MATRICES //
  //////////////////////


  inline size_t ColMajor::GetFirst(size_t i, size_t j)
  {
    return j;
  }
  inline size_t ColMajor::GetSecond(size_t i, size_t j)
  {
    return i;
  }
  inline size_t ColMajor::GetBeginLoop(size_t i)
  {
    return 0;
  }
  inline size_t ColMajor::GetEndLoop(size_t m, size_t n, size_t i)
  {
    return m;
  }


  inline size_t RowMajor::GetFirst(size_t i, size_t j)
  {
    return i;
  }
  inline size_t RowMajor::GetSecond(size_t i, size_t j)
  {
    return j;
  }
  inline size_t RowMajor::GetBeginLoop(size_t i)
  {
    return 0;
  }
  inline size_t RowMajor::GetEndLoop(size_t m, size_t n, size_t i)
  {
    return n;
  }


  /////////////
  // VECTORS //
  /////////////


  class VectFull
  {
  };


  class VectSparse
  {
  };


  class Collection
  {
  };


  class DenseSparseCollection
  {
  };


  class PETScSeq
  {
  };


  class PETScPar
  {
  };


  class PETScSeqDense
  {
  };


  class PETScMPIDense
  {
  };


  class PETScMPIAIJ
  {
  };


  ////////////
  // SPARSE //
  ////////////


  inline size_t ColSparse::GetFirst(size_t i, size_t j)
  {
    return j;
  }
  inline size_t ColSparse::GetSecond(size_t i, size_t j)
  {
    return i;
  }
  inline size_t ColSparse::GetBeginLoop(size_t i)
  {
    return 0;
  }
  inline size_t ColSparse::GetEndLoop(size_t m, size_t n, size_t i)
  {
    return m;
  }


  inline size_t RowSparse::GetFirst(size_t i, size_t j)
  {
    return i;
  }
  inline size_t RowSparse::GetSecond(size_t i, size_t j)
  {
    return j;
  }
  inline size_t RowSparse::GetBeginLoop(size_t i)
  {
    return 0;
  }
  inline size_t RowSparse::GetEndLoop(size_t m, size_t n, size_t i)
  {
    return n;
  }


  inline size_t ColSymSparse::GetFirst(size_t i, size_t j)
  {
    return j;
  }
  inline size_t ColSymSparse::GetSecond(size_t i, size_t j)
  {
    return i;
  }
  inline size_t ColSymSparse::GetBeginLoop(size_t i)
  {
    return i;
  }
  inline size_t ColSymSparse::GetEndLoop(size_t m, size_t n, size_t i)
  {
    return m;
  }


  inline size_t RowSymSparse::GetFirst(size_t i, size_t j)
  {
    return i;
  }
  inline size_t RowSymSparse::GetSecond(size_t i, size_t j)
  {
    return j;
  }
  inline size_t RowSymSparse::GetBeginLoop(size_t i)
  {
    return i;
  }
  inline size_t RowSymSparse::GetEndLoop(size_t m, size_t n, size_t i)
  {
    return n;
  }


  ///////////////
  // SYMMETRIC //
  ///////////////


  inline size_t ColSymPacked::GetFirst(size_t i, size_t j)
  {
    return j;
  }
  inline size_t ColSymPacked::GetSecond(size_t i, size_t j)
  {
    return i;
  }
  inline size_t ColSymPacked::GetBeginLoop(size_t i)
  {
    return i;
  }
  inline size_t ColSymPacked::GetEndLoop(size_t m, size_t n, size_t i)
  {
    return m;
  }


  inline size_t RowSymPacked::GetFirst(size_t i, size_t j)
  {
    return i;
  }
  inline size_t RowSymPacked::GetSecond(size_t i, size_t j)
  {
    return j;
  }
  inline size_t RowSymPacked::GetBeginLoop(size_t i)
  {
    return i;
  }
  inline size_t RowSymPacked::GetEndLoop(size_t m, size_t n, size_t i)
  {
    return n;
  }


  inline size_t ColSym::GetFirst(size_t i, size_t j)
  {
    return j;
  }
  inline size_t ColSym::GetSecond(size_t i, size_t j)
  {
    return i;
  }
  inline size_t ColSym::GetBeginLoop(size_t i)
  {
    return i;
  }
  inline size_t ColSym::GetEndLoop(size_t m, size_t n, size_t i)
  {
    return m;
  }


  inline size_t RowSym::GetFirst(size_t i, size_t j)
  {
    return i;
  }
  inline size_t RowSym::GetSecond(size_t i, size_t j)
  {
    return j;
  }
  inline size_t RowSym::GetBeginLoop(size_t i)
  {
    return i;
  }
  inline size_t RowSym::GetEndLoop(size_t m, size_t n, size_t i)
  {
    return n;
  }


  ///////////////
  // HERMITIAN //
  ///////////////


  inline size_t ColHerm::GetFirst(size_t i, size_t j)
  {
    return j;
  }
  inline size_t ColHerm::GetSecond(size_t i, size_t j)
  {
    return i;
  }
  inline size_t ColHerm::GetBeginLoop(size_t i)
  {
    return 0;
  }
  inline size_t ColHerm::GetEndLoop(size_t m, size_t n, size_t i)
  {
    return i+1;
  }


  inline size_t RowHerm::GetFirst(size_t i, size_t j)
  {
    return i;
  }
  inline size_t RowHerm::GetSecond(size_t i, size_t j)
  {
    return j;
  }
  inline size_t RowHerm::GetBeginLoop(size_t i)
  {
    return i;
  }
  inline size_t RowHerm::GetEndLoop(size_t m, size_t n, size_t i)
  {
    return n;
  }


  inline size_t ColHermPacked::GetFirst(size_t i, size_t j)
  {
    return j;
  }
  inline size_t ColHermPacked::GetSecond(size_t i, size_t j)
  {
    return i;
  }
  inline size_t ColHermPacked::GetBeginLoop(size_t i)
  {
    return 0;
  }
  inline size_t ColHermPacked::GetEndLoop(size_t m, size_t n, size_t i)
  {
    return i+1;
  }


  inline size_t RowHermPacked::GetFirst(size_t i, size_t j)
  {
    return i;
  }
  inline size_t RowHermPacked::GetSecond(size_t i, size_t j)
  {
    return j;
  }
  inline size_t RowHermPacked::GetBeginLoop(size_t i)
  {
    return i;
  }
  inline size_t RowHermPacked::GetEndLoop(size_t m, size_t n, size_t i)
  {
    return n;
  }


  ////////////////
  // TRIANGULAR //
  ////////////////


  inline size_t ColUpTriang::GetFirst(size_t i, size_t j)
  {
    return j;
  }
  inline size_t ColUpTriang::GetSecond(size_t i, size_t j)
  {
    return i;
  }
  inline size_t ColUpTriang::GetBeginLoop(size_t i)
  {
    return 0;
  }
  inline size_t ColUpTriang::GetEndLoop(size_t m, size_t n, size_t i)
  {
    return i+1;
  }
  inline bool ColUpTriang::UpLo()
  {
    return true;
  }


  inline size_t ColLoTriang::GetFirst(size_t i, size_t j)
  {
    return j;
  }
  inline size_t ColLoTriang::GetSecond(size_t i, size_t j)
  {
    return i;
  }
  inline size_t ColLoTriang::GetBeginLoop(size_t i)
  {
    return i;
  }
  inline size_t ColLoTriang::GetEndLoop(size_t m, size_t n, size_t i)
  {
    return m;
  }
  inline bool ColLoTriang::UpLo()
  {
    return false;
  }


  inline size_t RowUpTriang::GetFirst(size_t i, size_t j)
  {
    return i;
  }
  inline size_t RowUpTriang::GetSecond(size_t i, size_t j)
  {
    return j;
  }
  inline size_t RowUpTriang::GetBeginLoop(size_t i)
  {
    return i;
  }
  inline size_t RowUpTriang::GetEndLoop(size_t m, size_t n, size_t i)
  {
    return n;
  }
  inline bool RowUpTriang::UpLo()
  {
    return true;
  }


  inline size_t RowLoTriang::GetFirst(size_t i, size_t j)
  {
    return i;
  }
  inline size_t RowLoTriang::GetSecond(size_t i, size_t j)
  {
    return j;
  }
  inline size_t RowLoTriang::GetBeginLoop(size_t i)
  {
    return 0;
  }
  inline size_t RowLoTriang::GetEndLoop(size_t m, size_t n, size_t i)
  {
    return i+1;
  }
  inline bool RowLoTriang::UpLo()
  {
    return false;
  }


  inline size_t ColUpTriangPacked::GetFirst(size_t i, size_t j)
  {
    return j;
  }
  inline size_t ColUpTriangPacked::GetSecond(size_t i, size_t j)
  {
    return i;
  }
  inline size_t ColUpTriangPacked::GetBeginLoop(size_t i)
  {
    return 0;
  }
  inline size_t ColUpTriangPacked::GetEndLoop(size_t m, size_t n, size_t i)
  {
    return i+1;
  }
  inline bool ColUpTriangPacked::UpLo()
  {
    return true;
  }


  inline size_t ColLoTriangPacked::GetFirst(size_t i, size_t j)
  {
    return j;
  }
  inline size_t ColLoTriangPacked::GetSecond(size_t i, size_t j)
  {
    return i;
  }
  inline size_t ColLoTriangPacked::GetBeginLoop(size_t i)
  {
    return i;
  }
  inline size_t ColLoTriangPacked::GetEndLoop(size_t m, size_t n, size_t i)
  {
    return m;
  }
  inline bool ColLoTriangPacked::UpLo()
  {
    return false;
  }


  inline size_t RowUpTriangPacked::GetFirst(size_t i, size_t j)
  {
    return i;
  }
  inline size_t RowUpTriangPacked::GetSecond(size_t i, size_t j)
  {
    return j;
  }
  inline size_t RowUpTriangPacked::GetBeginLoop(size_t i)
  {
    return i;
  }
  inline size_t RowUpTriangPacked::GetEndLoop(size_t m, size_t n, size_t i)
  {
    return n;
  }
  inline bool RowUpTriangPacked::UpLo()
  {
    return true;
  }


  inline size_t RowLoTriangPacked::GetFirst(size_t i, size_t j)
  {
    return i;
  }
  inline size_t RowLoTriangPacked::GetSecond(size_t i, size_t j)
  {
    return j;
  }
  inline size_t RowLoTriangPacked::GetBeginLoop(size_t i)
  {
    return 0;
  }
  inline size_t RowLoTriangPacked::GetEndLoop(size_t m, size_t n, size_t i)
  {
    return i+1;
  }
  inline bool RowLoTriangPacked::UpLo()
  {
    return false;
  }


  inline size_t ColComplexSparse::GetFirst(size_t i, size_t j)
  {
    return j;
  }
  inline size_t ColComplexSparse::GetSecond(size_t i, size_t j)
  {
    return i;
  }
  inline size_t ColComplexSparse::GetBeginLoop(size_t i)
  {
    return 0;
  }


  inline size_t RowComplexSparse::GetFirst(size_t i, size_t j)
  {
    return i;
  }
  inline size_t RowComplexSparse::GetSecond(size_t i, size_t j)
  {
    return j;
  }
  inline size_t RowComplexSparse::GetBeginLoop(size_t i)
  {
    return 0;
  }


  inline size_t ColSymComplexSparse::GetFirst(size_t i, size_t j)
  {
    return j;
  }
  inline size_t ColSymComplexSparse::GetSecond(size_t i, size_t j)
  {
    return i;
  }
  inline size_t ColSymComplexSparse::GetBeginLoop(size_t i)
  {
    return i;
  }


  inline size_t RowSymComplexSparse::GetFirst(size_t i, size_t j)
  {
    return i;
  }
  inline size_t RowSymComplexSparse::GetSecond(size_t i, size_t j)
  {
    return j;
  }
  inline size_t RowSymComplexSparse::GetBeginLoop(size_t i)
  {
    return i;
  }

} // namespace Seldon.

#define SELDON_FILE_STORAGE_INLINE_CXX
#endif
