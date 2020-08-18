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
const TVar &TPZPetScMatrix<TVar>::GetVal( const int64_t row, const int64_t col ) const
{
#ifdef PZDEBUG
    if(row >=  this->Rows() || row<0 || col >=  this->Cols() || col<0) {
        std::cout << "TPZPetScMatrix::operator() "," Index out of bounds \n";
        DebugStop();
    }
#endif

    const PetscScalar *vals;
    const PetscInt *cols;
    PetscInt ncols = this->Cols();
    MatGetRow(fMat, row, &ncols, &cols, &vals);
    return *(&vals[col]);
}

template<class TVar>
inline int TPZPetScMatrix<TVar>::PutVal(const int64_t row,const int64_t col,const TVar & value )
{
    return MatSetValue(fMat, row, col, value, INSERT_VALUES);
}

template<class TVar>
PetscScalar &TPZPetScMatrix<TVar>::operator()( const int64_t row, const int64_t col)
{
#ifdef PZDEBUG
    if(row >=  this->Rows() || row<0 || col >=  this->Cols() || col<0) {
        std::cout << "TPZPetScMatrix<TVar>::operator() "," Index out of bounds\n";
        DebugStop();
    }
#endif
    const PetscScalar *vals;
    const PetscInt *cols;
    PetscInt ncols = this->Cols();
    MatGetRow(fMat, row, &ncols, &cols, &vals);
    auto val = const_cast< PetscScalar * > (vals);
    return *(&val[col]);
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

// MUMPS methods - Direct solver
template <class TVar>
int TPZPetScMatrix<TVar>::Decompose_Cholesky(std::list<int64_t> &singular) {
	if (  this->fDecomposed && this->fDecomposed != ECholesky) 
    {
        std::cout << "Decompose_Cholesky <Matrix already Decomposed>\n";
        DebugStop();
    }
    if (  this->fDecomposed ) return ECholesky;
	if (  this->Rows()!=this->Cols() ){
        std::cout << "Decompose_Cholesky <Matrix must be square>\n";
        DebugStop();
    } 
	
    // MPI function: Determines the size of the group associated with a communicator
    // PETSC_COMM_WORLD: communicator
    PetscMPIInt size;
    int ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);

    KSP ksp; // linear solver context
    PC pc; // preconditioner context

    ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr); //Creating the solver context
    ierr = KSPSetOperators(ksp,fMat,fMat);CHKERRQ(ierr); // Setting the matrices associated with the linear system

    ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr); // Getting the pc variable to modify the default options
    ierr = PCSetType(pc,PCCHOLESKY);CHKERRQ(ierr); // Setting that the decomposition method is LU
    PCFactorSetMatSolverType(pc,MATSOLVERMUMPS); // Setting the direct solver 
    PCFactorSetUpMatSolverType(pc);

    // Setting icntl MUMPS parameters - More details at the end of this file
    // Just an example:
    PCFactorGetMatrix(pc,&fMat);
    int icntl = 7; // ICNTL(7) computes a symmetric permutation in case of sequential analysis
    int ival = 2; // 2 : Approximate Minimum Fill (AMF) is used
    MatMumpsSetIcntl(fMat,icntl,ival);

    MatCholeskyFactor(fMat, NULL, NULL); // Decomposing matrix
    // However, the documentation says it's better to use the solver directly
    // ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);

	return ECholesky;
}

template <class TVar>
int TPZPetScMatrix<TVar>::Decompose_LU(std::list<int64_t> &singular)
{
	if (  this->fDecomposed && this->fDecomposed != ELU)
    {
        std::cout << "Decompose_LU <Matrix already Decomposed>\n";
        DebugStop();
    }
	if (  this->fDecomposed ) return ELU;
	if (  this->Rows()!=this->Cols() ){
        std::cout << "Decompose_LU <Matrix must be square>\n";
        DebugStop();
    } 
	
    // MPI function: Determines the size of the group associated with a communicator
    // PETSC_COMM_WORLD: communicator
    PetscMPIInt size;
    int ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);

    KSP ksp; // linear solver context
    PC pc; // preconditioner context

    ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr); //Creating the solver context
    ierr = KSPSetOperators(ksp,fMat,fMat);CHKERRQ(ierr); // Setting the matrices associated with the linear system

    ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr); // Getting the pc variable to modify the default options
    ierr = PCSetType(pc,PCLU);CHKERRQ(ierr); // Setting that the decomposition method is LU
    PCFactorSetMatSolverType(pc,MATSOLVERMUMPS); // Setting the direct solver 
    PCFactorSetUpMatSolverType(pc);

    // Setting icntl MUMPS parameters - More details at the end of this file
    // Just an example:
    PCFactorGetMatrix(pc,&fMat);
    int icntl = 7; // ICNTL(7) computes a symmetric permutation in case of sequential analysis
    int ival = 2; // 2 : Approximate Minimum Fill (AMF) is used
    MatMumpsSetIcntl(fMat,icntl,ival);

    MatLUFactor(fMat, NULL, NULL, NULL); // Decomposing matrix
    // However, the documentation says it's better to use the solver directly
    // ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);

	return ELU;
}

template class TPZPetScMatrix<STATE>;
template class TPZPetScMatrix<float>;

#endif

/*
    MUMPS ICNTL PARAMETERS:

• ICNTL(1) is the output stream for error messages
• ICNTL(2) is the output stream for diagnostic printing, statistics, and warning message
• ICNTL(3) is the output stream for global information, collected on the host
• ICNTL(4) is the level of printing for error, warning, and diagnostic messages
• ICNTL(5) controls the matrix input format
• ICNTL(6) permutes the matrix to a zero-free diagonal and/or scale the matrix
• ICNTL(7) computes a symmetric permutation in case of sequential analysis
• ICNTL(8) describes the scaling strategy
• ICNTL(9) computes the solution using A or AT
• ICNTL(10) applies the iterative refinement to the computed solution
• ICNTL(11) computes statistics related to an error analysis
• ICNTL(12) defines an ordering strategy for symmetric matrices
• ICNTL(13) controls the parallelism of the root node
• ICNTL(14) controls the percentage increase in the estimated working space
• ICNTL(16) controls the setting of the number of OpenMP threads
• ICNTL(18) defines the strategy for the distributed input matrix
• ICNTL(19) computes the Schur complement matrix
• ICNTL(20) determines the format (dense, sparse, or distributed) of the right-hand sides
• ICNTL(21) determines the distribution (centralized or distributed) of the solution vectors
• ICNTL(22) controls the in-core/out-of-core (OOC) factorization and solve
• ICNTL(23) corresponds to the maximum size of the working memory in MegaBytes that MUMPS can allocate per working processor
• ICNTL(24) controls the detection of “null pivot rows”
• ICNTL(25) allows the computation of a solution of a deficient matrix and also of a null space basis
• ICNTL(26) drives the solution phase if a Schur complement matrix
• ICNTL(27) controls the blocking size for multiple right-hand sides
• ICNTL(28) determines whether a sequential or parallel computation of the ordering is performed
• ICNTL(29) defines the parallel ordering tool to be used to compute the fill-in reducing permutation
• ICNTL(30) computes a user-specified set of entries in the inverse A−1 of the original matrix
• ICNTL(31) indicates which factors may be discarded during the factorization
• ICNTL(32) performs the forward elimination of the right-hand sides during the factorization
• ICNTL(33) computes the determinant of the input matrix
• ICNTL(34) controls the deletion of the files in case of save/restore
• ICNTL(35) controls the activation of the Block Low-Rank (BLR) feature
• ICNTL(36) controls the choice of BLR factorization variant
• ICNTL(37) reserved in current version
• ICNTL(38) estimated compression rate of LU factors
• ICNTL(39) reserved in current version
• ICNTL(40-57) not used in current version
• ICNTL(58) reserved in current version
• ICNTL(59-60) not used in current version
*/