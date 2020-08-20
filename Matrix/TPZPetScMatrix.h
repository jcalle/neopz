//
//  TPZPetScMatrix.hpp
//  PZ
//
//  Created by Karolinne Coelho on 11/08/2020.
//
//

#pragma once

#include "pz_config.h"

// #ifdef USING_PETSC

#define PETSC_USE_64BIT_INDICES

#include <stdio.h>
#include <mpi.h>
#include <petscksp.h>
#include "pzmanvector.h"
#include "tpzautopointer.h"
#include "pzmatrix.h"
#include "TPZStructMatrixBase.h"
#include "pzskylmat.h"

class TPZStructMatrixOR;

class TPZStructMatrixBase;

template<class TVar>
class TPZPetScMatrix: public TPZMatrix<TVar>, public TPZStructMatrix
{
public:
    enum MSystemType {ESparse, EDense, EBlock};
    
    enum MStructure {EStructureSymmetric, EStructureNonSymmetric};
    
    enum MProperty {EPositiveDefinite, EIndefinite};

    enum MMode {EParallel, ESerial};

    Mat fMat;

    inline TPZPetScMatrix() : TPZMatrix<TVar>(){
        fSystemType = ESparse;
        MatCreateSeqAIJ(PETSC_COMM_SELF, 0, 0, 0, NULL, &fMat);
        MatAssemblyBegin(this->fMat,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(this->fMat,MAT_FINAL_ASSEMBLY);
    }

    inline  TPZPetScMatrix (Mat m){
        fMat = m;
        MatGetSize(m, &(this->fRow), &(this->fCol));
        if((this->fRow)*(this->fCol)) fElem = new TVar[(this->fRow)*(this->fCol)];
        fSystemType = ESparse;
    }

    TPZPetScMatrix(TPZCompMesh * mesh) : TPZStructMatrix(mesh)
    {
    }
    
    TPZPetScMatrix(TPZAutoPointer<TPZCompMesh> cmesh) : TPZStructMatrix(cmesh)
    {
    }
	
    inline  TPZPetScMatrix (int rows, int columns) : TPZMatrix<TVar>(rows,columns){
        this->fRow = rows;
        this->fCol = columns;
        MatCreateSeqAIJ(PETSC_COMM_SELF, rows, columns, 0, NULL, &fMat);
        if(rows*columns) fElem = new TVar[rows*columns];
        MatAssemblyBegin(this->fMat,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(this->fMat,MAT_FINAL_ASSEMBLY);
    }

    inline  TPZPetScMatrix (int rows, int columns, MSystemType sys) :  TPZMatrix<TVar>(rows, columns){
        fSystemType = sys;
        if(rows*columns) fElem = new TVar[rows*columns];
        switch (sys)
        {
        case ESparse:
            MatCreateSeqAIJ(PETSC_COMM_SELF, rows, columns, 0, NULL, &fMat);
            MatAssemblyBegin(this->fMat,MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(this->fMat,MAT_FINAL_ASSEMBLY);
            break;

        case EDense:
            MatCreateSeqDense(PETSC_COMM_SELF, rows, columns, NULL, &fMat);
            MatAssemblyBegin(this->fMat,MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(this->fMat,MAT_FINAL_ASSEMBLY);
            break;
        
        case EBlock:
            std::cout <<"TPZPetScMatrix::EBlock - Type not implemented yet.\n";
            DebugStop();
            break;

        }
        this->fRow = rows;
        this->fCol = columns;
    }

    inline  TPZPetScMatrix (int64_t rows, int64_t columns, REAL val) : TPZMatrix<TVar>(rows,columns)
    {
        MatCreateSeqAIJ(PETSC_COMM_SELF, rows, columns, 0, NULL, &fMat);
        this->fRow = rows;
        this->fCol = columns;
        if(rows*columns) fElem = new TVar[rows*columns];
        const PetscInt idm = 0;
        const PetscInt idn = 0;
        MatSetOption(this->fMat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
        for (PetscInt i = 0; i < rows; i++)
        {
            for (PetscInt j = 0; j < columns; j++)
            {
                MatSetValues(this->fMat, 1, &i, 1, &j, &val, INSERT_VALUES);
            }
            
        }
        MatAssemblyBegin(this->fMat,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(this->fMat,MAT_FINAL_ASSEMBLY);
    }

    inline TPZPetScMatrix (int m, int n, MSystemType sys, MMode mode) : TPZPetScMatrix (m, n, sys)
    {
        if (mode == ESerial)
        {
            return;
        }
        
        fSystemType = sys;
        this->fRow = m;
        this->fCol = n;
        PetscMPIInt size;
        MPI_Comm_size(PETSC_COMM_WORLD,&size);

        // I'm not configuring this values, but according to the manual setting then can inprove the performance by more than 50x:
        PetscInt d_nz = m; //number of nonzeros per row in DIAGONAL portion of local submatrix (same value is used for all local rows)
        PetscInt o_nz = m; //number of nonzeros per row in the OFF-DIAGONAL portion of local submatrix (same value is used for all local rows).
        switch (sys)
        {
        case ESparse:
            MatCreateAIJ(PETSC_COMM_WORLD, m, n, PETSC_DETERMINE, PETSC_DETERMINE, d_nz, NULL, o_nz, NULL, &fMat);
            break;

        case EDense:
            MatCreateDense(PETSC_COMM_WORLD, m, n, PETSC_DECIDE, PETSC_DECIDE, NULL, &fMat);
            break;
        
        case EBlock:
            std::cout <<"TPZPetScMatrix:: Not implemented yet.\n";
            DebugStop();
            break;
        }
    }

    TPZPetScMatrix(const int64_t dim, const TPZVec<int64_t> &skyline ) :
    TPZRegisterClassId(&TPZPetScMatrix::ClassId),TPZMatrix<TVar>( dim, dim ), fElem(dim+1), fStorage(0)
    {
        fElem.Fill(0);
        InitializeElem(skyline,fStorage,fElem);
    }

    TPZPetScMatrix(const TPZPetScMatrix& cp)
       : fSystemType(cp.fSystemType), fMat(cp.fMat), fStructure(cp.fStructure), fProperty(cp.fProperty)
    {
        this->fRow = cp.Rows();
        this->fCol = cp.Cols();
        this->fMesh = cp.Mesh();
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
    }

    const TVar &GetVal(const int64_t row,const int64_t col ) const override;

    int PutVal(const int64_t row,const int64_t col,const TVar & value ) override;

    // Generic operators
    PetscScalar &operator()(const int64_t row,const int64_t col);

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

    // Direct solvers
    int Decompose_Cholesky();

    int Decompose_LU();

    int Subst_Forward( TPZFMatrix<TVar> *B ) const;

    int Subst_Backward( TPZFMatrix<TVar> *B ) const;
	
    virtual TPZMatrix<STATE> * Create()
    {
        TPZVec<int64_t> skyline;
        fMesh->Skyline(skyline);
        fEquationFilter.FilterSkyline(skyline);
        int64_t neq = fEquationFilter.NActiveEquations();
        return this->ReallyCreate(neq,skyline);//new TPZSkylMatrix<STATE>(neq,skyline);
    }

    TPZMatrix<STATE> * ReallyCreate(int64_t neq, const TPZVec<int64_t> &skyline){
        return new TPZSkylMatrix<STATE>(neq,skyline);
    }

    TPZStructMatrix * Clone(){
        return new TPZPetScMatrix<TVar>(*this);
    }

    void InitializeElem(const TPZVec<int64_t> &skyline, TPZVec<TVar> &storage, TPZVec<TVar *> &point)
    {
        int64_t dim = skyline.NElements();
        int64_t nel = NumElements(skyline);
        storage.Resize(nel);
        storage.Fill(0.);
        point.Resize(dim+1);
        if(dim) {
            point[0] = &storage[0];
            point[dim] = &storage[0]+nel;
        } 
        else {
            point[0] = 0;
        }
        
        for(int64_t i=1; i<dim+1; i++) { 
            point[i] = point[i-1] + (i-1)-skyline[i-1]+1;
        }
    }

    int64_t NumElements(const TPZVec<int64_t> &skyline)
{
    int64_t dim = skyline.NElements();
    
    int64_t nelem=0;
    for(int64_t i=0; i<dim; i++) {
        nelem += i-skyline[i]+1;
    }
    
    return nelem;
}


protected:
    
    MSystemType fSystemType;

    MStructure fStructure;

    MProperty fProperty;

    MMode fMode;

    MatType fMattype;

    TPZVec<TVar *> fElem;

    KSP fKsp; // linear solver context

	TPZVec<TVar> fStorage;
    
};

// #endif