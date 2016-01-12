// Copyright (C) 2001-2011 Vivien Mallet
// Copyright (C) 2001-2011 Marc Duruflé
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


// To be included by Seldon.hxx

#ifndef SELDON_FILE_MATRIX_SYMSPARSE_HXX

#include "../share/Common.hxx"
#include "../share/Properties.hxx"
#include "../share/Storage.hxx"
#include "../share/Errors.hxx"
#include "../share/Allocator.hxx"

namespace Seldon
{


  //! Symmetric sparse-matrix class.
  /*!
    Symmetric sparse matrices are defined by: (1) the number of rows
    and columns; (2) the number of non-zero entries; (3) an array 'ptr_'
    of start indices (i.e. indices of the first element of each row or column,
    depending on the storage); (4) an array 'ind_' of column or row indices
    of each non-zero entry; (5) values of non-zero entries.\par
    Only values of the upper part are stored.
  */
  template <class T, class Prop, class Storage, class Allocator
	    = typename SeldonDefaultAllocator<Storage, T>::allocator>
  class Matrix_SymSparse: public Matrix_Base<T, Allocator>
  {
    // typedef declaration.
  public:
    typedef typename Allocator::value_type value_type;
    typedef typename Allocator::pointer pointer;
    typedef typename Allocator::const_pointer const_pointer;
    typedef typename Allocator::reference reference;
    typedef typename Allocator::const_reference const_reference;
    typedef typename Allocator::value_type entry_type;
    typedef typename Allocator::value_type access_type;
    typedef typename Allocator::value_type const_access_type;
    typedef typename SeldonDefaultAllocator<VectFull, size_t>::allocator AllocatorInt;

    // Attributes.
  protected:
    // Number of non-zero (stored) elements.
    size_t nz_;
    // Index (in data_) of first element stored for each row or column.
    size_t* ptr_;
    // Column or row index (in the matrix) each element.
    size_t* ind_;

    // Methods.
  public:
    // Constructors.
    Matrix_SymSparse();
    Matrix_SymSparse(size_t i, size_t j);
    Matrix_SymSparse(size_t i, size_t j, size_t nz);
    template <class Storage0, class Allocator0,
	      class Storage1, class Allocator1,
	      class Storage2, class Allocator2>
    Matrix_SymSparse(size_t i, size_t j, Vector<T, Storage0, Allocator0>& values,
		     Vector<size_t, Storage1, Allocator1>& ptr,
		     Vector<size_t, Storage2, Allocator2>& ind);
    Matrix_SymSparse(const Matrix_SymSparse<T, Prop, Storage, Allocator>& A);

    // Destructor.
    ~Matrix_SymSparse();
    void Clear();

    // Memory management.
    template <class Storage0, class Allocator0,
	      class Storage1, class Allocator1,
	      class Storage2, class Allocator2>
    void SetData(size_t i, size_t j,
		 Vector<T, Storage0, Allocator0>& values,
		 Vector<size_t, Storage1, Allocator1>& ptr,
		 Vector<size_t, Storage2, Allocator2>& ind);
    void SetData(size_t i, size_t j, size_t nz, pointer values, size_t* ptr, size_t* ind);
    void Nullify();
    void Reallocate(size_t i, size_t j);
    void Reallocate(size_t i, size_t j, size_t nz);
    void Resize(size_t i, size_t j);
    void Resize(size_t i, size_t j, size_t nz);
    void Copy(const Matrix_SymSparse<T, Prop, Storage, Allocator>& A);

    // Basic methods.
    size_t GetNonZeros() const;
    size_t GetDataSize() const;
    int64_t GetMemorySize() const;
    size_t* GetPtr() const;
    size_t* GetInd() const;
    size_t GetPtrSize() const;
    size_t GetIndSize() const;

    // Element acess and affectation.
    const value_type operator() (size_t i, size_t j) const;
    value_type& Val(size_t i, size_t j);
    value_type& Get(size_t i, size_t j);
    const value_type& Val(size_t i, size_t j) const;
    const value_type& Get(size_t i, size_t j) const;
    void Set(size_t i, size_t j, const T& x);
    void AddInteraction(size_t i, size_t j, const T& x);

    void AddInteractionRow(size_t i, size_t nb, const Vector<size_t>& col,
			   const Vector<T>& val);

    Matrix_SymSparse<T, Prop, Storage, Allocator>&
    operator= (const Matrix_SymSparse<T, Prop, Storage, Allocator>& A);

    // Convenient functions.
    void Zero();
    void SetIdentity();
    void Fill();
    template <class T0>
    void Fill(const T0& x);
    void FillRand();
    
    void Print() const;
    void Write(string FileName) const;
    void Write(ostream& FileStream) const;
    void WriteText(string FileName, bool cplx = false) const;
    void WriteText(ostream& FileStream, bool cplx = false) const;
    void Read(string FileName);
    void Read(istream& FileStream);
    void ReadText(string FileName, bool cplx = false);
    void ReadText(istream& FileStream, bool cplx = false);

#ifdef SELDON_WITH_VIRTUAL
    typedef typename ClassComplexType<T>::Treal Treal;
    typedef typename ClassComplexType<T>::Tcplx Tcplx;

    virtual void ApplySor(Vector<T>& x, const Vector<T>& r,
			  const typename ClassComplexType<T>::Treal& omega,
			  int nb_iter, int stage_ssor) const;

    virtual void ApplySor(const class_SeldonTrans&, Vector<T>& x, const Vector<T>& r,
			  const typename ClassComplexType<T>::Treal& omega,
			  int nb_iter, int stage_ssor) const;

    virtual void MltAddVector(const Treal& alpha, const Vector<Treal>& x,
			      const Treal& beta, Vector<Treal>& y) const;

    virtual void MltAddVector(const Tcplx& alpha, const Vector<Tcplx>& x,
			      const Tcplx& beta, Vector<Tcplx>& y) const;

    virtual void MltAddVector(const Treal& alpha, const SeldonTranspose&,
			      const Vector<Treal>& x,
			      const Treal& beta, Vector<Treal>& y) const;

    virtual void MltAddVector(const Tcplx& alpha, const SeldonTranspose&,
			      const Vector<Tcplx>& x,
			      const Tcplx& beta, Vector<Tcplx>& y) const;
    
    virtual void MltVector(const Vector<Treal>& x, Vector<Treal>& y) const;
    virtual void MltVector(const Vector<Tcplx>& x, Vector<Tcplx>& y) const;
    
    virtual void MltVector(const SeldonTranspose&,
			   const Vector<Treal>& x, Vector<Treal>& y) const;

    virtual void MltVector(const SeldonTranspose&,
			   const Vector<Tcplx>& x, Vector<Tcplx>& y) const;
    
    virtual bool IsSymmetric() const;
#endif
    
  };


  //! Column-major sparse-matrix class.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ColSymSparse, Allocator>:
    public Matrix_SymSparse<T, Prop, ColSymSparse, Allocator>
  {
    // typedef declaration.
  public:
    typedef typename Allocator::value_type value_type;
    typedef Prop property;
    typedef ColSymSparse storage;
    typedef Allocator allocator;

  public:
    Matrix();
    explicit Matrix(size_t i, size_t j);
    explicit Matrix(size_t i, size_t j, size_t nz);
    template <class Storage0, class Allocator0,
	      class Storage1, class Allocator1,
	      class Storage2, class Allocator2>
    Matrix(size_t i, size_t j,
	   Vector<T, Storage0, Allocator0>& values,
	   Vector<size_t, Storage1, Allocator1>& ptr,
	   Vector<size_t, Storage2, Allocator2>& ind);
  };


  //! Row-major sparse-matrix class.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, RowSymSparse, Allocator>:
    public Matrix_SymSparse<T, Prop, RowSymSparse, Allocator>
  {
    // typedef declaration.
  public:
    typedef typename Allocator::value_type value_type;
    typedef Prop property;
    typedef RowSymSparse storage;
    typedef Allocator allocator;

  public:
    Matrix();
    explicit Matrix(size_t i, size_t j);
    explicit Matrix(size_t i, size_t j, size_t nz);
    template <class Storage0, class Allocator0,
	      class Storage1, class Allocator1,
	      class Storage2, class Allocator2>
    Matrix(size_t i, size_t j,
	   Vector<T, Storage0, Allocator0>& values,
	   Vector<size_t, Storage1, Allocator1>& ptr,
	   Vector<size_t, Storage2, Allocator2>& ind);
  };

} // namespace Seldon.

#define SELDON_FILE_MATRIX_SYMSPARSE_HXX
#endif
