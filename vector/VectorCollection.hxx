// Copyright (C) 2010, INRIA
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


#ifndef SELDON_FILE_VECTOR_VECTORCOLLECTION_HXX


#include <map>
#include "../share/Common.hxx"
#include "../share/Properties.hxx"
#include "../share/Storage.hxx"
#include "../share/Errors.hxx"
#include "../share/Allocator.hxx"


namespace Seldon
{


  //! Structure for distributed vectors.
  template <class T, class Allocator >
  class Vector<T, Collection, Allocator>: public Vector_Base<T, Allocator>
  {
    // typedef declarations.
  public:
    typedef typename T::value_type value_type;
    typedef typename T::pointer pointer;
    typedef typename T::const_pointer const_pointer;
    typedef typename T::reference reference;
    typedef typename T::const_reference const_reference;

    typedef T vector_type;
    typedef vector_type* vector_pointer;
    typedef const vector_type* const_vector_pointer;
    typedef vector_type& vector_reference;
    typedef const vector_type& const_vector_reference;

    typedef Vector<vector_type, VectFull, NewAlloc<vector_type> >
    collection_type;
    typedef const collection_type const_collection_type;
    typedef collection_type& collection_reference;
    typedef const collection_type& const_collection_reference;

    typedef Collection storage;

    // Attributes.
  protected:
    // Number of vectors.
    size_t Nvector_;
    // Lengths of the underlying vectors.
    Vector<size_t, VectFull, MallocAlloc<size_t> > length_;
    // Cumulative sum of the lengths of the underlying vectors.
    Vector<size_t, VectFull, MallocAlloc<size_t> > length_sum_;
    // Pointers of the underlying vectors.
    collection_type vector_;

    //! Indexes of the inner vectors that have a name.
    map<string, size_t> label_map_;
    //! Names associated with the inner vectors.
    vector<string> label_vector_;

    //! Methods.
  public:
    // Constructor.
    explicit Vector();
    explicit Vector(size_t i);
    Vector(const Vector<T, Collection, Allocator>& A);

    // Destructor.
    ~Vector();
    void Clear();
    void Reallocate(size_t i);
    void Deallocate();

   // Management of the vectors.
    template <class T0, class Storage0, class Allocator0>
    void AddVector(const Vector<T0, Storage0, Allocator0>& vector,
                   string name="", bool duplicate_data = false);
    template <class T0, class Storage0, class Allocator0>
    void SetVector(size_t i, const Vector<T0, Storage0, Allocator0>& vector,
                   string name = "", bool duplicate_data = false);
    template <class T0, class Storage0, class Allocator0>
    void SetVector(string name, const Vector<T0, Storage0,
                   Allocator0>& vector, bool duplicate_data = false);
    void SetName(size_t i, string name);

    void SetData(const Vector<T, Collection, Allocator>& X);
    void Nullify();

    // Basic methods.
    size_t GetM() const;
    size_t GetLength() const;
    size_t GetNvector() const;

    const Vector<size_t, VectFull, MallocAlloc<size_t> >& GetVectorLength() const;
    const Vector<size_t, VectFull, MallocAlloc<size_t> >& GetLengthSum() const;

    size_t GetVectorIndex(string name) const;
    size_t GetIndex(string name) const;

    collection_reference GetVector();
    const_collection_reference GetVector() const;
    vector_reference GetVector(size_t i);
    const_vector_reference GetVector(size_t i) const;
    vector_reference GetVector(string name);
    const_vector_reference GetVector(string name) const;

    // Element access and assignment.
    const_reference operator() (size_t i) const;
    reference operator() (size_t i);
    Vector<T, Collection, Allocator>& operator=
    (const Vector<T, Collection, Allocator>& X);

    void Copy(const Vector<T, Collection, Allocator>& X,
              bool duplicate_data = true);
    template <class T0, class Allocator0>
    void Copy(const Vector<T0, VectFull, Allocator0>& X,
              bool duplicate_data = true);

    template <class T0>
    Vector<T, Collection, Allocator>& operator*= (const T0& X);

    // Convenient method.
    template <class T0>
    void Fill(const T0& x);
    void Print() const;

    // Input/output functions.
    void Write(string FileName, bool with_size = true) const;
    void Write(ostream& FileStream, bool with_size = true) const;
    void WriteText(string FileName) const;
    void WriteText(ostream& FileStream) const;

    void Read(string FileName);
    void Read(string FileName, Vector<size_t, VectFull, MallocAlloc<size_t> >& length_);
    void Read(istream& FileStream, Vector<size_t, VectFull, MallocAlloc<size_t> >& length_);


  };


  template <class T, class Allocator>
  ostream& operator << (ostream& out,
			const Vector<T, Collection, Allocator>& V);


} // namespace Seldon.


#define SELDON_FILE_VECTOR_VECTORCOLLECTION_HXX
#endif
