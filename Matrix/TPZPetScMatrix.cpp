//
//  TPZPetScMatrix.hpp
//  PZ
//
//  Created by Karolinne Coelho on 11/08/2020.
//
//

#include "TPZPetScMatrix.h"

#ifdef USING_PETSC

#include <petscmat.h>

#include "pzsysmp.h"
#include "pzysmp.h"
#include "pzlog.h"
#include "pzmatrix.h"

template<class TVar>
inline const TVar &TPZPetScMatrix<TVar>::GetVal( const int64_t row, const int64_t col ) const
{
#ifdef PZDEBUG
    if(row >=  this->Rows() || row<0 || col >=  this->Cols() || col<0) {
        std::cout << "TPZPetScMatrix::operator() "," Index out of bounds \n";
        DebugStop();
    }
#endif
    PetscScalar value;
    MatGetValues(fMat, 1, &row, 1, &col, &value);
    return value;
}

template<class TVar>
inline int TPZPetScMatrix<TVar>::PutVal(const int64_t row,const int64_t col,const TVar & value )
{
    return MatSetValue(fMat, row, col, value, INSERT_VALUES);
}

template<class TVar>
inline TVar &TPZPetScMatrix<TVar>::operator()( const int64_t row, const int64_t col)
{
#ifndef NODEBUG
    if(row >=  this->Rows() || row<0 || col >=  this->Cols() || col<0) {
        std::cout << "TPZPetScMatrix<TVar>::operator() "," Index out of bounds\n";
        DebugStop();
    }
#endif
    // return *(&fMat+col*this->fRow+row); // don't know how to do it
}

// template<class TVar>
// &TPZPetScMatrix<TVar>::operator= (const TPZPetScMatrix<TVar> &A )
// {
//     fMat = A;
//     MatGetType(A, &fMattype);
//     MatGetSize(A, &(this->fRow), &(this->fCol));
// }



// template <class TVar>
// TPZPetScMatrix<TVar> TPZPetScMatrix<TVar>::operator+(const TPZPetScMatrix<TVar> &A ) const{
//     int ierr = MatSetValue(Mat m,PetscInt row,PetscInt col,PetscScalar value,InsertMode mode)
// }

template class TPZPetScMatrix<double>;
template class TPZPetScMatrix<long double>;
template class TPZPetScMatrix<float>;



#endif
