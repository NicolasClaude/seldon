%module seldon
%{
#include "SeldonHeader.hxx"
  %}

%include "std_string.i"

using namespace std;

// Include the header file with above prototypes
%include "SeldonHeader.hxx"
%include "Common/Common.hxx"
%include "Common/Storage.hxx"
%include "Common/Properties.hxx"
%include "Vector/Vector.hxx"
%include "Matrix/Matrix_Base.hxx"
%include "Matrix/Matrix_Pointers.hxx"
%include "Common/Allocator.hxx"

namespace Seldon
{
  %extend Vector<int, Vect_Full, MallocAlloc<int> >
  {
    int __getitem__(int index) {
      if (index < self->GetM())
	return self->GetData()[index];
      else
	return 0;
    }
    void __setitem__(int index, int value) {
      if (index >= 0 && index < self->GetM()) {
	self->GetData()[index] = value;
      }
    }
    unsigned long __len__() {
      return self->GetM();
    }
  }
  %extend Vector<double, Vect_Full, MallocAlloc<double> >
  {
    double __getitem__(int index) {
      if (index < self->GetM())
	return self->GetData()[index];
      else
	return 0;
    }
    void __setitem__(int index, double value) {
      if (index >= 0 && index < self->GetM()) {
	self->GetData()[index] = value;
      }
    }
    unsigned long __len__() {
      return self->GetM();
    }
  }

  %extend Matrix<int, General, RowMajor, MallocAlloc<int> >
  {
    int __getitem__(PyObject* args)
    {
      int i, j;
      int success = PyArg_ParseTuple(args, "ii", &i, &j);
      if (!success)
	throw std::out_of_range("Failed!");
      return (*self)(i, j);
    }
    Seldon::Vector<int, Seldon::Vect_Full, Seldon::MallocAlloc<int> > __getitem__(int i)
    {
      Seldon::Vector<int, Seldon::Vect_Full, Seldon::MallocAlloc<int> > v(self->GetN());
      for (int j = 0; j < self->GetN(); j++)
	v(j) = (*self)(i, j);
      return v;
    }
    void __setitem__(PyObject* args, int value)
    {
      int i, j;
      int success = PyArg_ParseTuple(args, "ii", &i, &j);
      if (!success)
	throw std::out_of_range("Failed!");
      (*self)(i, j) = value;
    }
    unsigned long __len__()
    {
      return self->GetM();
    }
  }
  %extend Matrix<double, General, RowMajor, MallocAlloc<double> >
  {
    double __getitem__(PyObject* args)
    {
      int i, j;
      int success = PyArg_ParseTuple(args, "ii", &i, &j);
      if (!success)
	throw std::out_of_range("Failed!");
      return (*self)(i, j);
    }
    Seldon::Vector<double, Seldon::Vect_Full, Seldon::MallocAlloc<double> > __getitem__(int i)
    {
      Seldon::Vector<double, Seldon::Vect_Full, Seldon::MallocAlloc<double> > v(self->GetN());
      for (int j = 0; j < self->GetN(); j++)
	v(j) = (*self)(i, j);
      return v;
    }
    void __setitem__(PyObject* args, double value)
    {
      int i, j;
      int success = PyArg_ParseTuple(args, "ii", &i, &j);
      if (!success)
	throw std::out_of_range("Failed!");
      (*self)(i, j) = value;
    }
    unsigned long __len__()
    {
      return self->GetM();
    }
  }

  %template(IntMalloc) MallocAlloc<int>;
  %template(BaseSeldonVectorInt) Vector_Base<int, MallocAlloc<int> >;
  %template(VectorInt) Vector<int, Vect_Full, MallocAlloc<int> >;
  %template(DoubleMalloc) MallocAlloc<double>;
  %template(BaseSeldonVectorDouble) Vector_Base<double, MallocAlloc<double> >;
  %template(VectorDouble) Vector<double, Vect_Full, MallocAlloc<double> >;

  %template(MatrixBaseInt) Matrix_Base<int, MallocAlloc<int> >;
  %template(MatrixPointersInt) Matrix_Pointers<int, General, RowMajor, MallocAlloc<int> >;
  %template(MatrixInt) Matrix<int, General, RowMajor, MallocAlloc<int> >;
  %template(MatrixBaseDouble) Matrix_Base<double, MallocAlloc<double> >;
  %template(MatrixPointersDouble) Matrix_Pointers<double, General, RowMajor, MallocAlloc<double> >;
  %template(MatrixDouble) Matrix<double, General, RowMajor, MallocAlloc<double> >;
}