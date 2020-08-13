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

#ifdef USING_PETSC

#define PETSC_USE_64BIT_INDICES

#include <stdio.h>
#include <petscksp.h>
#include "pzmanvector.h"
#include "tpzautopointer.h"
#include "pzmatrix.h"

template<class TVar>
class TPZSYsmpMatrix;

template<class TVar>
class TPZFYsmpMatrix;

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
            std::cout <<"TPZPetScMatrix:: Not implemented yet.\n";
            DebugStop();
            break;

        default:
            std::cout << "TPZPetScMatrix:: System type not defined.\n";
            DebugStop();
            break;

        }
    }

    inline TPZPetScMatrix (int m, int n, MSystemType sys, MMode mode) : TPZPetScMatrix (m, n, sys)
    {
        fSystemType = sys;
        this->fRow = m;
        this->fCol = n;
        MPI_Comm comm;
        switch (sys)
        {
        case ESparse:
            MatCreateSeqAIJ(comm, m, n, 0, NULL, &fMat);
            break;

        case EDense:
            MatCreateDense(comm, m, n, PETSC_DECIDE, PETSC_DECIDE, NULL, &fMat);
            break;
        
        case EBlock:
            std::cout <<"TPZPetScMatrix:: Not implemented yet.\n";
            DebugStop();
            break;

        default:
            std::cout << "TPZPetScMatrix:: System type not defined.\n";
            DebugStop();
            break;
        }
    }

    TPZPetScMatrix(TPZPetScMatrix& cp)
       : fSystemType(cp.fSystemType), fMat(cp.fMat), fStructure(cp.fStructure), fProperty(cp.fProperty)
    {
    }

    TPZPetScMatrix(TPZPetScMatrix&& mv)
        : fSystemType(mv.fSystemType), fMat(mv.fMat), fStructure(mv.fStructure), fProperty(mv.fProperty)
    {
    }

    const TVar &GetVal(const int64_t row,const int64_t col ) const override;

    int PutVal(const int64_t row,const int64_t col,const TVar & value ) override;

    // Generic operators
    TVar &operator()(const int64_t row,const int64_t col);
    
    // TVar &operator()(const int64_t row);

    virtual TPZPetScMatrix&operator= (const TPZPetScMatrix<TVar> &A );
    
    // virtual TPZPetScMatrix<TVar>& operator= (const std::initializer_list<TVar>& list);
	// TPZPetScMatrix<TVar>& operator= (const std::initializer_list< std::initializer_list<TVar> >& list);
    // TPZPetScMatrix<TVar> operator+  (const TPZPetScMatrix<TVar> &A ) const;
    // TPZPetScMatrix<TVar> operator-  (const TPZPetScMatrix<TVar> &A ) const;
    // TPZPetScMatrix<TVar> operator*  ( TPZPetScMatrix<TVar> A ) const ;
    // TPZPetScMatrix<TVar> &operator+=(const TPZPetScMatrix<TVar> &A );
    // TPZPetScMatrix<TVar> &operator-=(const TPZPetScMatrix<TVar> &A );

    // TPZPetScMatrix<TVar> &operator= (const TVar val );
    // TPZPetScMatrix<TVar> operator+  (const TVar val ) const;
    // TPZPetScMatrix<TVar> operator-  (const TVar val ) const;
    // TPZPetScMatrix<TVar> operator*  (const TVar val ) const;
    // TPZPetScMatrix<TVar> &operator+=(const TVar val );
    // TPZPetScMatrix<TVar> &operator-=(const TVar val )  { return operator+=( -val ); }
    // TPZPetScMatrix<TVar> &operator*=(const TVar val );

protected:
    
    MSystemType fSystemType;

    Mat fMat;

    MStructure fStructure;

    MProperty fProperty;

    MMode fMode;

    MatType fMattype;
    
};

#endif
#endif /* TPZPetScMatrix_hpp */
