/**
 * @file
 * @brief Contains the implementation of the TPZSkylineStructMatrix methods. 
 */

#include "pzskylstrmatrix.h"
#include "pzskylmat.h"
#include "pzcmesh.h"
#include "pzsubcmesh.h"
#include "pzvec.h"

TPZStructMatrix * TPZSkylineStructMatrix::Clone(){
    return new TPZSkylineStructMatrix(*this);
}

TPZSkylineStructMatrix::TPZSkylineStructMatrix(const TPZSkylineStructMatrix &cp)
:TPZStructMatrix(cp) {
	//nothing here
}

TPZSkylineStructMatrix::TPZSkylineStructMatrix(TPZCompMesh *mesh) : TPZStructMatrix(mesh)
{
}

TPZSkylineStructMatrix::TPZSkylineStructMatrix(TPZAutoPointer<TPZCompMesh> cmesh) : TPZStructMatrix(cmesh)
{
}

TPZMatrix<STATE> * TPZSkylineStructMatrix::CreateAssemble(TPZFMatrix<STATE> &rhs,TPZAutoPointer<TPZGuiInterface> guiInterface){
	int neq = fMesh->NEquations();
	if(HasRange()) neq = fMaxEq-fMinEq;
	TPZMatrix<STATE> *stiff = Create();
	rhs.Redim(neq,1);
	Assemble(*stiff,rhs, guiInterface);
    return stiff;
}

TPZMatrix<STATE> * TPZSkylineStructMatrix::Create(){
    int neq = fMesh->NEquations();
    TPZVec<int> skyline;
	if (fOnlyInternal) {
		TPZSubCompMesh *submesh = dynamic_cast<TPZSubCompMesh *> (fMesh);
		if (submesh) {
			submesh->SkylineInternal(skyline);
		}
		else {
			fMesh->Skyline(skyline);
		}
	}
	else {
		fMesh->Skyline(skyline);
	}
    if(HasRange())
    {
		neq = fMaxEq-fMinEq;
		FilterSkyline(skyline);
    }
	else {
		// This statement is needed for compatibility with TPZSubCompMesh The number of equations of the stiffness matrix corresponds to "only" the internal nodes
		neq = skyline.NElements();
	}
	
    return this->ReallyCreate(neq,skyline);//new TPZSkylMatrix<STATE>(neq,skyline);
}


TPZMatrix<STATE> * TPZSkylineStructMatrix::ReallyCreate(int neq, const TPZVec<int> &skyline){
    return new TPZSkylMatrix<STATE>(neq,skyline);
}




TPZSkylineStructMatrix::~TPZSkylineStructMatrix(){
}

/*!
 \fn TPZSkylineStructMatrix::FilterSkyline()
 */
void TPZSkylineStructMatrix::FilterSkyline(TPZVec<int> &skyline)
{
	if(!HasRange()) return;
	//  int neq = skyline.NElements();
	int ieq;
	for(ieq = fMinEq; ieq < fMaxEq; ieq++)
	{
		skyline[ieq-fMinEq] = skyline[ieq]-fMinEq;
	}
	skyline.Resize(fMaxEq-fMinEq);
}
