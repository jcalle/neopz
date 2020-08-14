//
//  TPZPetScMatrix.hpp
//  PZ
//
//  Created by Karolinne Coelho on 11/08/2020.
//
//

#ifndef TPZPetScMatrix_hpp
#define TPZPetScMatrix_hpp

#include "pz_config.h"

// #ifdef USING_PETSC

#define PETSC_USE_64BIT_INDICES

#include <stdio.h>
#include <mpi.h>
#include <petscksp.h>
#include "pzmanvector.h"
#include "tpzautopointer.h"
#include "pzmatrix.h"

template<class TVar>
class TPZPetScMatrix: public TPZMatrix<TVar>
{
public:
    enum MSystemType {ESparse, EDense, EBlock};
    
    enum MStructure {EStructureSymmetric, EStructureNonSymmetric};
    
    enum MProperty {EPositiveDefinite, EIndefinite};

    enum MMode {EParallel, ESerial};

    inline TPZPetScMatrix(){
        fSystemType = ESparse;
        this->fRow = 0;
        this->fCol = 0;
        MatCreateSeqAIJ(PETSC_COMM_SELF, 0, 0, 0, NULL, &fMat);
    }

    inline  TPZPetScMatrix (Mat m){
        fMat = m;
        MatGetSize(m, &(this->fRow), &(this->fCol));
        fSystemType = ESparse;
    }

    inline  TPZPetScMatrix (int m, int n){
        this->fRow = m;
        this->fCol = n;
        MatCreateSeqAIJ(PETSC_COMM_SELF, m, n, 0, NULL, &fMat);
    }

    inline  TPZPetScMatrix (int m, int n, MSystemType sys){
        fSystemType = sys;
        switch (sys)
        {
        case ESparse:
            MatCreateSeqAIJ(PETSC_COMM_SELF, m, n, 0, NULL, &fMat);
            break;

        case EDense:
            MatCreateSeqDense(PETSC_COMM_SELF, m, n, NULL, &fMat);
            break;
        
        case EBlock:
            std::cout <<"TPZPetScMatrix::EBlock - Type not implemented yet.\n";
            DebugStop();
            break;

        default:
            std::cout << "TPZPetScMatrix:: System type not defined.\n";
            DebugStop();
            break;

        }
        this->fRow = m;
        this->fCol = n;
    }

    inline  TPZPetScMatrix (int64_t m, int64_t n, REAL val)
    {
        MatCreateSeqAIJ(PETSC_COMM_SELF, m, n, 0, NULL, &fMat);
        this->fRow = m;
        this->fCol = n;
        const PetscInt idm = int(m)-1;
        const PetscInt idn = int(n)-1;
        MatSetOption(this->fMat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
        MatSetValues(this->fMat, m, &idm, n, &idn, &val, INSERT_VALUES);
        MatAssemblyBegin(this->fMat,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(this->fMat,MAT_FINAL_ASSEMBLY);
    }

    // inline TPZPetScMatrix (int m, int n, MSystemType sys, MMode mode) : TPZPetScMatrix (m, n, sys)
    // {
    //     fSystemType = sys;
    //     this->fRow = m;
    //     this->fCol = n;
    //     switch (sys)
    //     {
    //     case ESparse:
    //         MatCreateSeqAIJ(PETSC_COMM_WORLD, m, n, 0, NULL, &fMat);
    //         break;

    //     case EDense:
    //         MatCreateDense(PETSC_COMM_WORLD, m, n, PETSC_DECIDE, PETSC_DECIDE, NULL, &fMat);
    //         break;
        
    //     case EBlock:
    //         std::cout <<"TPZPetScMatrix:: Not implemented yet.\n";
    //         DebugStop();
    //         break;

    //     default:
    //         std::cout << "TPZPetScMatrix:: System type not defined.\n";
    //         DebugStop();
    //         break;
    //     }
    // }

    TPZPetScMatrix(const TPZPetScMatrix& cp)
       : fSystemType(cp.fSystemType), fMat(cp.fMat), fStructure(cp.fStructure), fProperty(cp.fProperty)
    {
        this->fRow = cp.Rows();
        this->fCol = cp.Cols();
    }

    TPZPetScMatrix(TPZPetScMatrix&& mv)
        : fSystemType(mv.fSystemType), fMat(mv.fMat), fStructure(mv.fStructure), fProperty(mv.fProperty)
    {
        this->fRow = std::move(mv.fRow);
        this->fCol = std::move(mv.fCol);
    }

    CLONEDEF(TPZPetScMatrix<TVar>)
    TPZPetScMatrix(const TPZMatrix<TVar> & cp);

    virtual  ~TPZPetScMatrix()
    {
        // MatDestroy(&(this->fMat));
    }

    const TVar &GetVal(const int64_t row,const int64_t col ) const override;

    int PutVal(const int64_t row,const int64_t col,const TVar & value ) override;

    // Generic operators
    TVar &operator()(const int64_t row,const int64_t col);

    virtual TPZPetScMatrix<TVar>& operator= (const TPZPetScMatrix<TVar> &A );
    TPZPetScMatrix<TVar> operator+  (const TPZPetScMatrix<TVar> &A ) const;
    TPZPetScMatrix<TVar> operator-  (const TPZPetScMatrix<TVar> &A ) const;
    TPZPetScMatrix<TVar> operator*  ( TPZPetScMatrix<TVar> A ) const ;
    TPZPetScMatrix<TVar> &operator+=(const TPZPetScMatrix<TVar> &A );
    TPZPetScMatrix<TVar> &operator-=(const TPZPetScMatrix<TVar> &A );

    TPZPetScMatrix<TVar> operator+  (const TVar val ) const;
    TPZPetScMatrix<TVar> operator-  (const TVar val ) const
    {
        return operator+(-val);
    };
    TPZPetScMatrix<TVar> operator* (const TVar val ) const;
    TPZPetScMatrix<TVar> &operator+=(const TVar val );
    TPZPetScMatrix<TVar> &operator-=(const TVar val ) 
    { 
        return operator+=( -val );
    }
    TPZPetScMatrix<TVar> &operator*=(const TVar val );

    /** @brief Redimension a matrix, but maintain your elements. */
    int Resize(const int64_t newRows,const int64_t wCols ) override;
    
    /** @brief Redimension a matrix and ZERO your elements. */
    int Redim(const int64_t newRows,const int64_t newCols ) override;
    
    /** @brief Makes Zero all the elements */
    int Zero() override;

protected:
    
    MSystemType fSystemType;

    Mat fMat;

    MStructure fStructure;

    MProperty fProperty;

    MMode fMode;

    MatType fMattype;
    
};

// #endif
#endif /* TPZPetScMatrix_hpp */
