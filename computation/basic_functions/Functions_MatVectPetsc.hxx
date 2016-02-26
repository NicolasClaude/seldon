#ifndef SELDON_FILE_FUNCTIONS_MATVECT_PETSC_HXX
#define SELDON_FILE_FUNCTIONS_MATVECT_PETSC_HXX


namespace Seldon
{


  ////////////
  /// MLT  ///


  template <class T1, class Prop1, class Allocator1,
	    class T2, class Allocator2,
	    class T4, class Allocator4>
  void MltVector(const Matrix<T1, Prop1, PETScMPIDense, Allocator1>& M,
		 const Vector<T2, PETScPar, Allocator2>& X,
		 Vector<T4, PETScPar, Allocator4>& Y);


  /// MLT  ///
  ////////////


  ////////////
  // MLTADD //


  template <class T0,
            class T1, class Prop1, class Allocator1,
            class T2, class Allocator2,
            class T3,
            class T4, class Allocator4>
  void MltAddVector(const T0 alpha,
              const Matrix<T1, Prop1, PETScMPIAIJ, Allocator1>& M,
              const Vector<T2, PETScSeq, Allocator2>& X,
              const T3 beta, Vector<T4, PETScSeq, Allocator4>& Y);

  template <class T0,
            class T1, class Prop1, class Storage1, class Allocator1,
            class T2, class Storage2, class Allocator2,
            class T3,
            class T4, class Allocator4>
  void MltAddVector(const T0 alpha,
              const Matrix<T1, Prop1, Storage1, Allocator1>& M,
              const Vector<T2, Storage2, Allocator2>& X,
              const T3 beta, Vector<T4, PETScPar, Allocator4>& Y);

  template <class T0,
            class Prop1, class Allocator1,
            class Allocator2,
            class Allocator4>
  void MltAddVector(const T alpha,
              const Matrix<T, Prop1, PETScMPIAIJ, Allocator1>& M,
              const Vector<T, PETScPar, Allocator2>& X,
              const T beta, Vector<T, PETScPar, Allocator4>& Y);

  template <class T0,
            class T1, class Prop1, class Allocator1,
            class T2, class Allocator2,
            class T3,
            class T4, class Allocator4>
  void MltAddVector(const T0 alpha,
              const Matrix<T1, Prop1, PETScMPIDense, Allocator1>& M,
              const Vector<T2, PETScPar, Allocator2>& X,
              const T3 beta, Vector<T4, PETScPar, Allocator4>& Y);

  template <class T0,
            class T1, class Prop1, class Allocator1,
            class T2, class Allocator2,
            class T3,
            class T4, class Allocator4>
  void MltAddVector(const T0 alpha,
              const Matrix<T1, Prop1, PETScMPIAIJ, Allocator1>& M,
              const Vector<T2, VectFull, Allocator2>& X,
              const T3 beta, Vector<T4, PETScSeq, Allocator4>& Y);

  template <class T0,
            class T1, class Prop1, class Allocator1,
            class T2, class Allocator2,
            class T3,
            class T4, class Allocator4>
  void MltAddVector(const T0 alpha,
              const Matrix<T1, Prop1, PETScMPIAIJ, Allocator1>& M,
              const Vector<T2, VectFull, Allocator2>& X,
              const T3 beta, Vector<T4, PETScPar, Allocator4>& Y);

  template <class T0,
            class T1, class Prop1, class Allocator1,
            class T2, class Allocator2,
            class T3,
            class T4, class Allocator4>
  void MltAddVector(const T0 alpha,
              const Matrix<T1, Prop1, PETScMPIDense, Allocator1>& M,
              const Vector<T2, VectFull, Allocator2>& X,
              const T3 beta, Vector<T4, PETScPar, Allocator4>& Y);

  template <class T0,
            class T1, class Prop1, class Allocator1,
            class T2, class Allocator2,
            class T3,
            class T4, class Allocator4>
  void MltAddVector(const T0 alpha,
              const Matrix<T1, Prop1, PETScMPIAIJ, Allocator1>& M,
              const Vector<T2, PETScPar, Allocator2>& X,
              const T3 beta, Vector<T4, VectFull, Allocator4>& Y);

  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Allocator2,
	    class T3,
	    class T4, class Allocator4>
  void MltAddVector(const T0 alpha,
	      const SeldonTranspose& Trans,
	      const Matrix<T1, Prop1, PETScMPIAIJ, Allocator1>& M,
	      const Vector<T2, PETScPar, Allocator2>& X,
	      const T3 beta,
	      Vector<T4, PETScPar, Allocator4>& Y);

  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Allocator2,
	    class T3,
	    class T4, class Allocator4>
  void MltAddVector(const T0 alpha,
	      const SeldonTranspose& Trans,
	      const Matrix<T1, Prop1, PETScMPIDense, Allocator1>& M,
	      const Vector<T2, PETScPar, Allocator2>& X,
	      const T3 beta,
	      Vector<T4, PETScPar, Allocator4>& Y);

  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Allocator2,
	    class T3,
	    class T4, class Allocator4>
  void MltAddVector(const T0 alpha,
	      const SeldonTranspose& Trans,
	      const Matrix<T1, Prop1, PETScMPIAIJ, Allocator1>& M,
	      const Vector<T2, PETScPar, Allocator2>& X,
	      const T3 beta,
	      Vector<T4, VectFull, Allocator4>& Y);


  // MLTADD //
  ////////////


  /////////////
  // SOLVELU //


  template <class T0, class Prop0, class Allocator0,
	    class T1, class Allocator1>
  void GetAndSolveLU(Matrix<T0, Prop0, PETScMPIDense, Allocator0>& M,
                     Vector<T1, PETScPar, Allocator1>& Y);

  template <class T0, class Prop0,  class Allocator0,
	    class T1, class Allocator1>
  void SolveLuVector(const Matrix<T0, Prop0, PETScMPIDense, Allocator0>& M,
		     Vector<T1, PETScPar, Allocator1>& Y);


  // SOLVELU //
  /////////////

} // namespace Seldon.


#endif
