#ifndef SELDON_FILE_SELDONPETSC_H
#define SELDON_FILE_SELDONPETSC_H



#define IFMAINPROC                                                      \
  int rank_macro;                                                       \
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank_macro); if (rank_macro == 0) {

#define MPIDISP(x)  MPIIF( DISP(x) )

#define MPIERR(x)  MPIIF( ERR(x) )

#define MPIIF(x) IFMAINPROC x ;}

#include "vector/PetscVector.cxx"
#include "matrix/PetscMatrix.cxx"

#include "matrix/FunctionsPetsc.cxx"
#include "computation/basic_functions/Functions_MatrixPetsc.cxx"
#include "computation/basic_functions/Functions_MatVectPetsc.cxx"
#include "computation/basic_functions/Functions_VectorPetsc.cxx"

#endif
