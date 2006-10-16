/* Generated by Together */

#ifndef TPZBLOCKDIAGONALSTRUCTMATRIX_H
#define TPZBLOCKDIAGONALSTRUCTMATRIX_H
#include "pzstrmatrix.h"
class TPZCompMesh;
class TPZFMatrix;
class TPZMatrix;
class TPZBlockDiagonal;
template <class T>
class TPZVec;

/**
 * Implements Block Diagonal Structural Matrices
 * @ingroup structural
 */
class TPZBlockDiagonalStructMatrix : public TPZStructMatrix {
public:    

enum MBlockStructure {ENodeBased, EVertexBased, EElementBased};

  TPZBlockDiagonalStructMatrix(TPZCompMesh *);
  
  // create a sparse blockdiagonal matrix, overlapping should be assumed
  virtual TPZMatrix * Create();
  
  // create a sparse blockdiagonal matrix of the given color
  // this should be used to create a sequence solver
  TPZMatrix * Create(int color);
  
  virtual TPZMatrix * CreateAssemble(TPZFMatrix &rhs);

  virtual TPZStructMatrix * Clone();    

public:

    void AssembleBlockDiagonal(TPZBlockDiagonal & block);
private:
    void BlockSizes(TPZVec < int > & blocksizes);
    
    MBlockStructure fBlockStructure;
    int fOverlap;
    
    
};
#endif //TPZBLOCKDIAGONALSTRUCTMATRIX_H
