/**
 * @file
 * @brief Contains the TPZPetscSSpStructMatrix class which implements sparse structural matrices.
 */

#pragma once

#include "pzstrmatrix.h"
#include "pzysmp.h"

#include "pzcmesh.h"
#include "pzsubcmesh.h"
#include "pzelmat.h"

#include "pzmatrix.h"
#include "pzfmatrix.h"

/**
 * @brief Implements Sparse Structural Matrices. \ref structural "Structural Matrix"
 * @ingroup structural
 */
class TPZPetscSSpStructMatrix : public TPZStructMatrix {
    
public:    
	
    TPZPetscSSpStructMatrix(TPZCompMesh * mesh) : TPZStructMatrix(mesh)
    {
    }

    
    virtual ~TPZPetscSSpStructMatrix() override
    {
    }
	
    virtual TPZMatrix<STATE> * Create() override;
	
    virtual TPZMatrix<STATE> * SetupMatrixData(TPZStack<int64_t> & elgraph, TPZVec<int64_t> &elgraphindex);
    
    using TPZStructMatrix::CreateAssemble;
	virtual TPZMatrix<STATE> * CreateAssemble(TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) override;
	
    virtual TPZStructMatrix * Clone() override;
	
private :
    
    TPZPetscSSpStructMatrix(): TPZStructMatrix(){
    
    };
    
    friend TPZPersistenceManager;
};