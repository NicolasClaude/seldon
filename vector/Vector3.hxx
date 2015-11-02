// Copyright (C) 2010-2012, INRIA
// Author(s): Marc Fragu, Vivien Mallet
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


#ifndef SELDON_FILE_VECTOR_VECTOR_3_HXX


#ifndef SELDON_VECTOR3_DEFAULT_ALLOCATOR_0
/*! \def SELDON_VECTOR3_DEFAULT_ALLOCATOR_0
  Default allocator for the inner vectors in Vector3.
*/
#define SELDON_VECTOR3_DEFAULT_ALLOCATOR_0 SELDON_DEFAULT_ALLOCATOR
#endif

#ifndef SELDON_VECTOR3_DEFAULT_ALLOCATOR_1
/*! \def SELDON_VECTOR3_DEFAULT_ALLOCATOR_1
  Default allocator for the inner vector of vectors in Vector3.
*/
#define SELDON_VECTOR3_DEFAULT_ALLOCATOR_1 MallocObject
#endif

#ifndef SELDON_VECTOR3_DEFAULT_ALLOCATOR_2
#define SELDON_VECTOR3_DEFAULT_ALLOCATOR_2 MallocObject
/*! \def SELDON_VECTOR3_DEFAULT_ALLOCATOR_2
  Default allocator for the vector of vectors of vectors in Vector3.
*/
#endif


namespace Seldon
{

  //! %Vector of vectors of vectors.
  /*! Vector3 is a structure that acts like a vector of vectors of
    vectors. Both inner vectors and inner vectors of vectors can be of any
    dimension, so that this structure is more flexible than an Array3D.
    \tparam T numerical type of the inner vectors.
    \tparam Allocator0 allocator for the inner vectors. The default allocator
    is SELDON_DEFAULT_ALLOCATOR.
    \tparam Allocator1 allocator for the vector of vectors. It is recommended
    to choose NewAlloc or, for more efficient in reallocations, MallocObject
    (default allocator here): these allocators can manage an array of inner
    vectors.
    \tparam Allocator2 allocator for the vector of vectors of vectors. It is
    recommended to choose NewAlloc or, for more efficient in reallocations,
    MallocObject (default allocator here).
  */
  template <class T,
	    class Allocator0 = SELDON_VECTOR3_DEFAULT_ALLOCATOR_0<T>,
	    class Allocator1 = SELDON_VECTOR3_DEFAULT_ALLOCATOR_1<
	      Vector<T, Vect_Full, Allocator0> >,
	    class Allocator2 =
	    SELDON_VECTOR3_DEFAULT_ALLOCATOR_2<
	      Vector<Vector<T, Vect_Full, Allocator0>,
		     Vect_Full, Allocator1> > >
  class Vector3
  {
  public:
    typedef T value_type;
    typedef T* pointer;
    typedef const T* const_pointer;
    typedef T& reference;
    typedef const T& const_reference;

  protected:
    Vector<Vector<Vector<T, Vect_Full, Allocator0>, Vect_Full, Allocator1>,
	   Vect_Full, Allocator2> data_;

  public:

    /*** Constructors and destructor ***/

    Vector3();
    Vector3(size_t);
    Vector3(Vector<size_t>& length);
    template <class Allocator>
    Vector3(Vector<Vector<size_t>, Vect_Full, Allocator>& length);
    ~Vector3();

    /*** Management of the vectors ***/

    size_t GetLength() const;
    size_t GetSize() const;
    size_t GetLength(size_t i) const;
    size_t GetSize(size_t i) const;
    size_t GetLength(size_t i, size_t j) const;
    size_t GetSize(size_t i, size_t j) const;
    int64_t GetMemorySize() const;
    size_t GetNelement() const;
    size_t GetNelement(size_t beg, size_t end) const;
    size_t GetNelement(size_t beg0, size_t end0, size_t beg1, size_t end1) const;
    Vector<size_t> GetShape(size_t i) const;
    void GetShape(size_t i, Vector<size_t>& shape) const;
    void Reallocate(size_t N);
    void Reallocate(size_t i, size_t N);
    void Reallocate(size_t i, size_t j, size_t N);

    template <class Td, class Allocatord>
    void Flatten(Vector<Td, VectFull, Allocatord>& data) const;
    template <class Td, class Allocatord>
    void Flatten(size_t beg, size_t end, Vector<Td, VectFull, Allocatord>& data)
      const;
    template <class Td, class Allocatord>
    void Flatten(size_t beg0, size_t end0, size_t beg1, size_t end1,
                 Vector<Td, VectFull, Allocatord>& data) const;

    void PushBack(size_t i, size_t j, const T& x);
    void PushBack(size_t i, const Vector<T, Vect_Full, Allocator0>& X);
    void PushBack(const Vector<Vector<T, Vect_Full, Allocator0>,
		  Vect_Full, Allocator1>& X);
    void PushBack(const Vector<Vector<Vector<T, Vect_Full, Allocator0>,
                  Vect_Full, Allocator1>, Vect_Full, Allocator2>& X);
    void PushBack(const Vector3<T, Allocator0, Allocator1, Allocator2>& X);


    void Clear();
    void Clear(size_t i);
    void Clear(size_t i, size_t j);

    void Fill(const T& x);

    Vector<Vector<Vector<T, Vect_Full, Allocator0>, Vect_Full, Allocator1>,
	   Vect_Full, Allocator2>&
    GetVector();
    const Vector<Vector<Vector<T, Vect_Full, Allocator0>,
                        Vect_Full, Allocator1>, Vect_Full, Allocator2>&
    GetVector() const;

    Vector<Vector<T, Vect_Full, Allocator0>, VectFull, Allocator1>&
    GetVector(size_t i);
    const Vector<Vector<T, Vect_Full, Allocator0>, VectFull, Allocator1>&
    GetVector(size_t i) const;

    Vector<T, Vect_Full, Allocator0>& GetVector(size_t i, size_t j);
    const Vector<T, Vect_Full, Allocator0>& GetVector(size_t i, size_t j) const;

    /*** Element access and assignment ***/
    const
    Vector<Vector<T, Vect_Full, Allocator0>, VectFull, Allocator1>&
    operator() (size_t i) const;
    Vector<Vector<T, Vect_Full, Allocator0>, VectFull, Allocator1>&
    operator() (size_t i);

    const Vector<T, Vect_Full, Allocator0>& operator() (size_t i, size_t j)
      const;
    Vector<T, Vect_Full, Allocator0>& operator() (size_t i, size_t j);

    const_reference operator() (size_t i, size_t j, size_t k) const;
    reference operator() (size_t i, size_t j, size_t k);

    /*** Convenient method ***/

    void Print() const;

    /*** Input/output functions ***/

    void Write(string file_name, bool with_size = true) const;
    void Write(ostream& file_stream, bool with_size = true) const;
    void Read(string file_name, bool with_size = true);
    void Read(istream& file_stream, bool with_size = true);
  };

}


#define SELDON_FILE_VECTOR_VECTOR_3_HXX
#endif
