//
//  TPZPetScMatrix.hpp
//  PZ
//
//  Created by Karolinne Coelho on 11/08/2020.
//
//

#include "TPZPetScMatrix.h"

// #ifdef USING_PETSC

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
#ifdef PZDEBUG
    if(row >=  this->Rows() || row<0 || col >=  this->Cols() || col<0) {
        std::cout << "TPZPetScMatrix<TVar>::operator() "," Index out of bounds\n";
        DebugStop();
    }
#endif
    // return *(&fMat+col*this->fRow+row); // don't know how to do it
}

template<class TVar>
TPZPetScMatrix<TVar> &TPZPetScMatrix<TVar>::operator= (const TPZPetScMatrix<TVar> &A )
{   
#ifdef PZDEBUG
    if ((this->Rows() != A.Rows()) || (this->Cols() != A.Cols()) ) {
		std::cout << "Add(TPZMatrix<>&, TPZMatrix) <different dimensions>\n";
	}
#endif
    fMat = A.fMat;
    MatGetType(A.fMat, &fMattype);
}

template <class TVar>
TPZPetScMatrix<TVar> TPZPetScMatrix<TVar>::operator+(const TPZPetScMatrix<TVar> &A ) const
{
    if ( (A.Rows() != this->Rows())  ||  (A.Cols() != this->Cols()) )
    {
        std::cout << "Operator+ <matrices with different dimensions>\n";
    }
    TPZPetScMatrix<TVar> B(A);
    MatAXPY(B.fMat, 1, this->fMat, DIFFERENT_NONZERO_PATTERN); 

    return B;
}

template <class TVar>
TPZPetScMatrix<TVar> TPZPetScMatrix<TVar>::operator-(const TPZPetScMatrix<TVar> &A ) const
{
    if ( (A.Rows() != this->Rows())  ||  (A.Cols() != this->Cols()) )
    {
        std::cout << "Operator- <matrices with different dimensions>\n";
    }
    TPZPetScMatrix<TVar> B(A);
    MatAXPY(B.fMat, -1, this->fMat, DIFFERENT_NONZERO_PATTERN); 
    
    return B;
}

template <class TVar>
TPZPetScMatrix<TVar> TPZPetScMatrix<TVar>::operator*(TPZPetScMatrix<TVar> A ) const
{
    if ( this->Cols() != A.Rows() )
    {
        std::cout << "Operator- <matrices with different dimensions>\n";
    }
    MatMatMult(A.fMat, this->fMat, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &(A.fMat));
    // MatGetSize(this->fMat, &(this->fRow), &(this->fCol));

    return A;
}

template <class TVar>
TPZPetScMatrix<TVar> &TPZPetScMatrix<TVar>::operator+=(const TPZPetScMatrix<TVar> &A )
{
    if ( (A.Rows() != this->Rows())  ||  (A.Cols() != this->Cols()) )
    {
        std::cout << "Operator+ <matrices with different dimensions>\n";
    }
    MatAXPY(this->fMat, 1, A.fMat, DIFFERENT_NONZERO_PATTERN); 

    return *this;
}

template <class TVar>
TPZPetScMatrix<TVar> &TPZPetScMatrix<TVar>::operator-=(const TPZPetScMatrix<TVar> &A )
{
    if ( (A.Rows() != this->Rows())  ||  (A.Cols() != this->Cols()) )
    {
        std::cout << "Operator- <matrices with different dimensions>\n";
    }
    MatAXPY(this->fMat, -1, A.fMat, DIFFERENT_NONZERO_PATTERN); 
    
    return *this;
}

template <class TVar>
TPZPetScMatrix<TVar> TPZPetScMatrix<TVar>::operator+ (const TVar val) const
{
    TPZPetScMatrix<TVar> A(*this);
    auto vald = double(val);
    MatSetValues(A.fMat, A.fRow, &(A.fRow), A.fCol, &(A.fCol), &vald, ADD_VALUES);

    return A;
}

template <class TVar>
TPZPetScMatrix<TVar> TPZPetScMatrix<TVar>::operator* (const TVar val) const
{
    TPZPetScMatrix<TVar> A(*this);
    MatScale(A.fMat, val);

    return A;
}

template <class TVar>
TPZPetScMatrix<TVar> &TPZPetScMatrix<TVar>::operator+= (const TVar val)
{
    auto vald = double(val);
    MatSetValues(this->fMat, this->fRow, &(this->fRow), this->fCol, &(this->fCol), &vald, ADD_VALUES);
    return *this;
}

template <class TVar>
int TPZPetScMatrix<TVar>::Redim(const int64_t newRows,const int64_t newCols )
{
    int ierr = MatSetSizes(this->fMat, PETSC_DECIDE, PETSC_DECIDE, newRows, newCols);
    CHKERRQ(ierr);
    this->Zero();

    return ierr;
}

template <class TVar>
int TPZPetScMatrix<TVar>::Resize(const int64_t newRows,const int64_t newCols )
{
    int ierr = MatSetSizes(this->fMat, PETSC_DECIDE, PETSC_DECIDE, newRows, newCols);
    CHKERRQ(ierr);

    return ierr;
}

template <class TVar>
int TPZPetScMatrix<TVar>::Zero()
{
    MatZeroEntries(this->fMat);
}


template class TPZPetScMatrix<STATE>;
template class TPZPetScMatrix<float>;

// #endif