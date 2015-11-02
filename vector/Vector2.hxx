// Copyright (C) 2010, INRIA
// Author(s): Marc Fragu, Vivien Mallet
// Copyright (C) 2011, Vivien Mallet
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


#ifndef SELDON_FILE_VECTOR_VECTOR2_HXX


#ifndef SELDON_VECTOR2_DEFAULT_ALLOCATOR_0
/*! \def SELDON_VECTOR2_DEFAULT_ALLOCATOR_0
  Default allocator for the inner vectors in Vector2.
*/
#define SELDON_VECTOR2_DEFAULT_ALLOCATOR_0 SELDON_DEFAULT_ALLOCATOR
#endif

#ifndef SELDON_VECTOR2_DEFAULT_ALLOCATOR_1
/*! \def SELDON_VECTOR2_DEFAULT_ALLOCATOR_1
  Default allocator for the vector of vectors in Vector2.
*/
#define SELDON_VECTOR2_DEFAULT_ALLOCATOR_1 MallocObject
#endif


namespace Seldon
{


  //! %Vector of vectors.
  /*! Vector2 is a structure that acts like a vector of vectors. The inner
    vectors can be of any dimension, so that this structure is more flexible
    than a matrix.
    \tparam T numerical type of the inner vectors.
    \tparam Allocator0 allocator for the inner vectors. The default allocator
    is SELDON_DEFAULT_ALLOCATOR.
    \tparam Allocator1 allocator for the vector of vectors. It is recommended
    to choose NewAlloc or, for more efficient in reallocations, MallocObject
    (default allocator here): these allocators can manage an array of inner
    vectors.
  */
  template <class T,
            class Allocator0 = SELDON_VECTOR2_DEFAULT_ALLOCATOR_0<T>,
            class Allocator1 = SELDON_VECTOR2_DEFAULT_ALLOCATOR_1<
              Vector<T, VectFull, Allocator0> > >
  class Vector2
  {
  public:
    typedef T value_type;
    typedef T* pointer;
    typedef const T* const_pointer;
    typedef T& reference;
    typedef const T& const_reference;

  protected:
    Vector<Vector<T, VectFull, Allocator0>, VectFull, Allocator1> data_;

  public:

    /*** Constructors and destructor ***/

    Vector2();
    Vector2(size_t length);
    Vector2(const Vector<size_t>& length);
    ~Vector2();

    /*** Management of the vectors ***/

    bool IsEmpty() const;
    size_t GetLength() const;
    size_t GetSize() const;
    size_t GetLength(size_t i) const;
    size_t GetSize(size_t i) const;
    int64_t GetMemorySize() const;
    size_t GetNelement() const;
    size_t GetNelement(size_t beg, size_t end) const;
    Vector<size_t> GetShape() const;
    void GetShape(Vector<size_t>& shape) const;
    void Reallocate(size_t M);
    void Reallocate(size_t i, size_t N);
    void Reallocate(const Vector<size_t>& length);
    void Select(size_t beg, size_t end);
    Vector<T, VectFull, Allocator0> Flatten() const;
    template <class Td, class Allocatord>
    void Flatten(Vector<Td, VectFull, Allocatord>& data) const;
    template <class Td, class Allocatord>
    void Flatten(size_t beg, size_t end, Vector<Td, VectFull, Allocatord>& data)
      const;

    void PushBack(size_t i, const T& x);
    void PushBack(const Vector<T, VectFull, Allocator0>& X);
#ifndef SWIG
    void PushBack(const Vector<Vector<T, VectFull, Allocator0>,
		  VectFull, Allocator1>& V);
#endif
    void PushBack(const Vector2<T, Allocator0, Allocator1>& V);

    void Clear();
    void Clear(size_t i);

    void Fill(const T& x);

    Vector<Vector<T, VectFull, Allocator0>, VectFull, Allocator1>&
    GetVector();
#ifndef SWIG
    const Vector<Vector<T, VectFull, Allocator0>, VectFull,
                 Allocator1> GetVector() const;
#endif

    Vector<T, VectFull, Allocator0>& GetVector(size_t i);
#ifndef SWIG
    const Vector<T, VectFull, Allocator0>& GetVector(size_t i) const;
#endif

    void Copy(const Vector2<T, Allocator0, Allocator1>& V);
    Vector2<T, Allocator0, Allocator1> Copy() const;

    /*** Element access and assignment ***/

#ifndef SWIG
    const Vector<T, VectFull, Allocator0>& operator() (size_t i) const;
#endif
    Vector<T, VectFull, Allocator0>& operator() (size_t i);
#ifndef SWIG
    const_reference operator() (size_t i, size_t j) const;
#endif
    reference operator() (size_t i, size_t j);

    /*** Convenient methods ***/

    template <class V2>
    bool HasSameShape(const V2& V) const;
    void Print() const;

    /*** Input/output functions ***/

    void Write(string file_name, bool with_size = true) const;
    void Write(ostream& file_stream, bool with_size = true) const;
    void Read(string file_name, bool with_size = true);
    void Read(istream& file_stream, bool with_size = true);
  };


}


#define SELDON_FILE_VECTOR_VECTOR2_HXX
#endif
