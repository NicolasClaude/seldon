// Copyright (C) 2003-2011 Marc Duruflé
// Copyright (C) 2001-2011 Vivien Mallet
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


#ifndef SELDON_FILE_SPARSE_VECTOR_INLINE_CXX

#include "SparseVector.hxx"

namespace Seldon
{


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    On exit, the vector is empty.
  */
  template <class T, class Allocator>
  inline Vector<T, VectSparse, Allocator>::Vector():
    Vector<T, VectFull, Allocator>()
  {
    index_ = NULL;
  }


  //! Main constructor.
  /*! Builds a vector of a given size.
    \param i length of the vector.
  */
  template <class T, class Allocator>
  inline Vector<T, VectSparse, Allocator>::Vector(int i):
    Vector<T, VectFull, Allocator>(i)
  {

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	this->index_ = index_allocator_.allocate(i, this);

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->index_ = NULL;
	this->data_ = NULL;
      }

    if (this->index_ == NULL)
      {
	this->m_ = 0;
	this->data_ = NULL;
      }

    if (this->data_ == NULL && i != 0)
      throw NoMemory("Vector<VectSparse>::Vector(int)",
		     string("Unable to allocate memory for a vector of size ")
		     + to_str(i * sizeof(T)) + " bytes ("
		     + to_str(i) + " elements).");
#endif

  }


  //! Copy constructor.
  /*! Builds a copy of a vector.
    \param V vector to be copied.
  */
  template <class T, class Allocator>
  inline Vector<T, VectSparse, Allocator>::
  Vector(const Vector<T, VectSparse, Allocator>& V) :
    Vector<T, VectFull, Allocator>()
  {
    this->index_ = NULL;
    Copy(V);
  }


  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor.
  template <class T, class Allocator>
  inline Vector<T, VectSparse, Allocator>::~Vector()
  {
    Clear();
  }
  

  /*********************
   * MEMORY MANAGEMENT *
   *********************/


  //! Clears the vector.
  /*!
    Destructs the vector.
    \warning On exit, the vector is an empty vector.
  */
  template <class T, class Allocator>
  inline void Vector<T, VectSparse, Allocator>::Clear()
  {
    // 'data_' is released.
#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif
	if (this->data_ != NULL)
	  {
	    this->vect_allocator_.deallocate(this->data_, this->m_);
	    this->data_ = NULL;
	  }

	if (index_ != NULL)
	  {
	    index_allocator_.deallocate(index_, this->m_);
	    index_ = NULL;
	  }

	this->m_ = 0;

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->data_ = NULL;
	index_ = NULL;
	this->m_ = 0;
	return;
      }
#endif

  }


  //! Vector reallocation.
  /*!
    The vector is resized.
    \param i new length of the vector.
    \warning Depending on your allocator, previous non-zero entries may be
    lost.
  */
  template <class T, class Allocator>
  inline void Vector<T, VectSparse, Allocator>::Reallocate(int i)
  {

    if (i != this->m_)
      {

	this->m_ = i;

#ifdef SELDON_CHECK_MEMORY
	try
	  {
#endif

	    this->data_ =
	      reinterpret_cast<pointer>(this->vect_allocator_
					.reallocate(this->data_, i, this));

	    index_
	      = reinterpret_cast<int*>(this->index_allocator_
				       .reallocate(index_, i, this));

#ifdef SELDON_CHECK_MEMORY
	  }
	catch (...)
	  {
	    this->m_ = 0;
	    this->data_ = NULL;
	    this->index_ = NULL;
	    return;
	  }
	if (this->data_ == NULL)
	  {
	    this->m_ = 0;
	    this->index_ = NULL;
	    return;
	  }
#endif

      }
  }


  //! Changes the number of non-zero entries of the vector.
  /*! Changes the number of non-zero entries to \a n. If \a n non-zero entries
    are available before resizing, they are all kept. Otherwise, only the
    first \n non-zero entries are kept.
    \param n new number of non-zero entries of the vector.
  */
  template <class T, class Allocator>
  inline void Vector<T, VectSparse, Allocator>::Resize(int n)
  {

    if (n == this->m_)
      return;

    Vector<T, VectFull, Allocator> new_value(n);
    Vector<int> new_index(n);
    int Nmin = min(this->m_, n);
    for (int i = 0; i < Nmin; i++)
      {
	new_value(i) = this->data_[i];
	new_index(i) = index_[i];
      }

    SetData(new_value, new_index);
  }


  /*! \brief Changes the length of the vector and sets its data array (low
    level method). */
  /*!
    Reallocates a vector and sets the new data array. It is useful to create
    a vector from pre-existing data.
    \param i new length of the vector.
    \param data the new data array. \a data contains the new elements of the
    vector and must therefore contain \a i elements.
    \param index the new index array. \a index contains the new indices of the
    non-zero entries and it must therefore contain \a i elements.
    \warning \a data has to be used carefully outside the object.  Unless you
    use 'Nullify', \a data will be freed by the destructor, which means that
    \a data must have been allocated carefully. The vector allocator should be
    compatible.
    \note This method should only be used by advanced users.
  */
  template <class T, class Allocator>
  inline void Vector<T, VectSparse, Allocator>
  ::SetData(int i, T* data, int* index)
  {
    this->Clear();

    this->m_ = i;

    this->data_ = data;
    this->index_ = index;
  }


  /*! \brief Changes the length of the vector and sets its data array (low
    level method). */
  /*!
    Reallocates a vector and sets the new data array. It is useful to create
    a vector from pre-existing data.
    \param data the new data array. \a data contains the new elements of the
    vector and must therefore contain \a i elements.
    \param index the new index array. \a index contains the new indices of the
    non-zero entries and it must therefore contain \a i elements.
    \note Vectors \a data and \a index are empty vector on exit.
  */
  template <class T, class Allocator>
  template<class Allocator2>
  inline void Vector<T, VectSparse, Allocator>
  ::SetData(Vector<T, VectFull, Allocator2>& data, Vector<int>& index)
  {

#ifdef SELDON_CHECK_BOUNDS
    if (data.GetM() != index.GetM())
      throw WrongDim("Vector<VectSparse>::SetData ",
		     string("The data vector and the index vector should")
		     + " have the same size.\n  Size of the data vector: "
		     + to_str(data.GetM()) + "\n  Size of index vector: "
		     + to_str(index.GetM()));
#endif

    SetData(data.GetM(), data.GetData(), index.GetData());
    data.Nullify();
    index.Nullify();
  }


  /*! \brief Lets the current vector point to the data of a second vector (low
    level method). */
  /*! Reallocates the current vector and lets its data point to those of \a V.
    \param V the vector to which the current vector points to (on exit).
    \warning On exit, both \a V and the current vector point to the same
    arrays in memory. Only one of them should eventually deallocate the memory
    blocks. The other one should be nullified by the user. In case the current
    vector is responsible for the deallocations, its allocator should be
    compatible with the allocator that created the memory blocks (which is
    probably the allocator of \a V).
  */
  template <class T, class Allocator>
  template<class Allocator2>
  inline void Vector<T, VectSparse, Allocator>
  ::SetData(const Vector<T, VectSparse, Allocator2>& V)
  {
    SetData(V.GetM(), V.GetData(), V.GetIndex());
  }


  //! Clears the vector without releasing memory.
  /*!
    On exit, the vector is empty and the memory has not been released.
    It is useful for low level manipulations on a Vector instance.
    \warning Memory is not released.
  */
  template <class T, class Allocator>
  inline void Vector<T, VectSparse, Allocator>::Nullify()
  {
    this->m_ = 0;
    this->data_ = NULL;
    this->index_ = NULL;
  }


  /**********************************
   * ELEMENT ACCESS AND AFFECTATION *
   **********************************/


  //! Access operator.
  /*!
    \param i index.
    \return The value of the non-zero element #\a i.
  */
  template <class T, class Allocator>
  inline typename Vector<T, VectSparse, Allocator>::reference
  Vector<T, VectSparse, Allocator>::Value(int i)
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongIndex("Vector<VectSparse>::Value(int)",
		       string("Index should be in [0, ") + to_str(this->m_-1)
		       + "], but is equal to " + to_str(i) + ".");
#endif

    return this->data_[i];
  }


  //! Access operator.
  /*!
    \param i index.
    \return The value of the non-zero element #\a i.
  */
  template <class T, class Allocator>
  inline typename Vector<T, VectSparse, Allocator>::const_reference
  Vector<T, VectSparse, Allocator>::Value(int i) const
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongIndex("Vector<VectSparse>::Value(int)",
		       string("Index should be in [0, ") + to_str(this->m_-1)
		       + "], but is equal to " + to_str(i) + ".");
#endif

    return this->data_[i];
  }


  //! Access operator.
  /*!
    \param i index.
    \return The index of the non-zero element #\a i.
  */
  template <class T, class Allocator>
  inline int& Vector<T, VectSparse, Allocator>::Index(int i)
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongIndex("Vector<VectSparse>::Index(int)",
		       string("Index should be in [0, ") + to_str(this->m_-1)
		       + "], but is equal to " + to_str(i) + ".");
#endif

    return this->index_[i];
  }


  //! Access operator.
  /*!
    \param i index.
    \return The row number of the non-zero element #\a i.
  */
  template <class T, class Allocator>
  inline int Vector<T, VectSparse, Allocator>::Index(int i) const
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongIndex("Vector<VectSparse>::Index(int)",
		       string("Index should be in [0, ") + to_str(this->m_-1)
		       + "], but is equal to " + to_str(i) + ".");
#endif

    return this->index_[i];
  }


  //! Access operator.
  /*!
    \param i index.
    \return The value of the vector at \a i.
  */
  template <class T, class Allocator>
  inline typename Vector<T, VectSparse, Allocator>::value_type
  Vector<T, VectSparse, Allocator>::operator() (int i) const
  {
    int k = 0;
    T zero;
    SetComplexZero(zero);
    // Searching for the entry.
    while (k < this->m_ && index_[k] < i)
      k++;

    if (k >= this->m_ || index_[k] != i)
      // The entry does not exist, a zero is returned.
      return zero;

    return this->data_[k];
  }


  //! Access method.
  /*! Returns the value of element \a i.
    \param[in] i index.
    \return Element \a i of the vector.
    \throw WrongArgument No reference can be returned because the element is a
    zero entry (not stored in the vector).
  */
  template <class T, class Allocator>
  inline typename Vector<T, VectSparse, Allocator>::value_type
  Vector<T, VectSparse, Allocator>::operator() (int i)
  {
    int k = 0;
    T zero;
    SetComplexZero(zero);
    // Searching for the entry.
    while (k < this->m_ && index_[k] < i)
      k++;

    if (k >= this->m_ || index_[k] != i)
      // The entry does not exist, a zero is returned.
      return zero;
    
    return this->data_[k];
  }
  
  
  //! Access method.
  /*! Returns the value of element \a i.
    \param[in] i index.
    \return Element \a i of the vector.
    \throw WrongArgument No reference can be returned because the element is a
    zero entry (not stored in the vector).
  */
  template <class T, class Allocator>
  inline typename Vector<T, VectSparse, Allocator>::reference
  Vector<T, VectSparse, Allocator>::Get(int i)
  {
    int k = 0;
    T zero;
    SetComplexZero(zero);
    // Searching for the entry.
    while (k < this->m_ && index_[k] < i)
      k++;

    if (k >= this->m_ || index_[k] != i)
      // The entry does not exist yet, so a zero entry is introduced.
      AddInteraction(i, zero);
    
    return this->data_[k];
  }


  //! Access method.
  /*! Returns the value of element \a i.
    \param[in] i index.
    \return Element \a i of the vector.
    \throw WrongArgument No reference can be returned because the element is a
    zero entry (not stored in the vector).
  */
  template <class T, class Allocator>
  inline typename Vector<T, VectSparse, Allocator>::const_reference
  Vector<T, VectSparse, Allocator>::Get(int i) const
  {
    int k = 0;
    // Searching for the entry.
    while (k < this->m_ && index_[k] < i)
      k++;

    if (k >= this->m_ || index_[k] != i)
      throw WrongArgument("Vector<VectSparse>::Val(int)",
                          "No reference to element " + to_str(i)
                          + " can be returned: it is a zero entry.");
    
    return this->data_[k];
  }


  //! Access method.
  /*! Returns the value of element \a i.
    \param[in] i index.
    \return Element \a i of the vector.
    \throw WrongArgument if i does not belong to the sparsity pattern
  */
  template <class T, class Allocator>
  inline typename Vector<T, VectSparse, Allocator>::reference
  Vector<T, VectSparse, Allocator>::Val(int i)
  {
    int k = 0;
    // Searching for the entry.
    while (k < this->m_ && index_[k] < i)
      k++;

    if (k >= this->m_ || index_[k] != i)
      throw WrongArgument("Vector<VectSparse>::Val(int)",
                          "the entry " + to_str(i) +
                          " does not belong to the sparsity pattern.");

    
    return this->data_[k];
  }


  //! Access method.
  /*! Returns the value of element \a i.
    \param[in] i index.
    \return Element \a i of the vector.
    \throw WrongArgument if i does not belong to the sparsity pattern
  */
  template <class T, class Allocator>
  inline typename Vector<T, VectSparse, Allocator>::const_reference
  Vector<T, VectSparse, Allocator>::Val(int i) const
  {
    int k = 0;
    // Searching for the entry.
    while (k < this->m_ && index_[k] < i)
      k++;

    if (k >= this->m_ || index_[k] != i)
      throw WrongArgument("Vector<VectSparse>::Val(int)",
                          "the entry " + to_str(i) +
                          " does not belong to the sparsity pattern.");
    
    return this->data_[k];
  }

  
  //! Duplicates a vector (assignment operator).
  /*!
    \param X vector to be copied.
    \note Memory is duplicated: \a X is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Allocator>
  inline Vector<T, VectSparse, Allocator>& Vector<T, VectSparse, Allocator>
  ::operator= (const Vector<T, VectSparse, Allocator>& X)
  {
    this->Copy(X);

    return *this;
  }


  //! Duplicates a vector.
  /*!
    \param X vector to be copied.
    \note Memory is duplicated: \a X is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Allocator>
  inline void Vector<T, VectSparse, Allocator>
  ::Copy(const Vector<T, VectSparse, Allocator>& X)
  {
    this->Reallocate(X.GetLength());

    this->vect_allocator_.memorycpy(this->data_, X.GetData(), this->m_);
    this->index_allocator_.memorycpy(this->index_, X.GetIndex(), this->m_);
  }


  /*******************
   * BASIC FUNCTIONS *
   *******************/


  /*! \brief Returns a pointer to the array containing the indices of the
    non-zero entries. */
  /*!
    \return A pointer to the array of the indices of the non-zero entries.
  */
  template <class T, class Allocator>
  inline int* Vector<T, VectSparse, Allocator>::GetIndex() const
  {
    return this->index_;
  }

  
} // namespace Seldon.

#define SELDON_FILE_SPARSE_VECTOR_INLINE_CXX
#endif