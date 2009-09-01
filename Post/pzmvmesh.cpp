#include "pzmvmesh.h"
#include "pzcmesh.h"
#include "pzmaterial.h"
#include "pzgraphnode.h"
#include "pzgraphel.h"
using namespace std;

TPZMVGraphMesh::TPZMVGraphMesh(TPZCompMesh *cmesh, int dimension, TPZAutoPointer<TPZMaterial> mat) : TPZGraphMesh(cmesh, dimension, mat) {
	fNumCases = 0;
	fNumSteps = 0;
   fStyle = EMVStyle;
}

TPZMVGraphMesh::TPZMVGraphMesh(TPZCompMesh *cmesh, int dimension, TPZMVGraphMesh *graph,TPZAutoPointer<TPZMaterial> mat) :
                TPZGraphMesh(cmesh, dimension,mat) {
  if(!mat) fMaterial = graph->fMaterial;
	fNumCases = graph->fNumCases;
	fNumSteps = graph->fNumSteps;
  fStyle = EMVStyle;
}

void TPZMVGraphMesh::DrawMesh(int numcases) {

	fNumCases = numcases;
	fNumSteps = 0;
	(fOutFile) << "%HEADER" << endl;
	(fOutFile) << "Arquivo gerado por PZ" << endl;
	(fOutFile) << "%RESULT\n1" << endl;
	(fOutFile) << 1 << " 'Caso Unico'" << endl;
	(fOutFile) << "%RESULT.CASE" << endl;
	(fOutFile) << 1 << " " << numcases << endl;
	for(int i=0; i<numcases;i++) {
		(fOutFile) << (i+1) << " 'step" << i << "'" << endl;
	}
	DrawNodes();
	DrawConnectivity(ECube);
}

void TPZMVGraphMesh::DrawSolution(int step, REAL time){

	fNumSteps++;
	int numscal = this->fScalarNames.NElements();
	int numvec = fVecNames.NElements();
	TPZVec<int> scalind(0);
	TPZVec<int> vecind(0);
	scalind.Resize(numscal);
	vecind.Resize(numvec);
	scalind.Fill(-1,0,numscal);
	vecind.Fill(-1,0,numvec);
	TPZAutoPointer<TPZMaterial> matp = Material();
   if(!matp) {
   	cout << "TPZMVGraphMesh no material found\n";
   	return;
	}
	int n;
	for(n=0; n<numscal; n++) {
		scalind[n] = matp->VariableIndex( fScalarNames[n]);
	}
	for(n=0; n<numvec; n++) {
		vecind[n] = matp->VariableIndex(fVecNames[n]);
	}

	(fOutFile) << "%RESULT.CASE.STEP\n" << (step+1) << endl;
	(fOutFile) << "%RESULT.CASE.STEP.TIME\n" << time << endl;
	int dispind = matp->VariableIndex("Displacement6");
	TPZVec<REAL> disp(6,0.);

	(fOutFile) << "%RESULT.CASE.STEP.NODAL.DISPLACEMENT" << endl;
	(fOutFile) << NPoints() << " 'Nodal Displ'" << endl;
   int nnod = fNodeMap.NElements(),i;
   for(i=0;i<nnod;i++) {
		TPZGraphNode *n = &fNodeMap[i];
		if(n) n->DrawSolution(dispind, EMVStyle);
	}
	(fOutFile) << "%RESULT.CASE.STEP.NODAL.SCALAR" << endl;
	(fOutFile) << numscal << endl;
	for(n=0;n<numscal;n++) {
		(fOutFile) << "'" <<   fScalarNames[n] << "' ";
	}
	(fOutFile) << endl;

	(fOutFile) << "%RESULT.CASE.STEP.NODAL.SCALAR.DATA" << endl;
	(fOutFile) << NPoints() << endl;
   nnod = fNodeMap.NElements();
   for(i=0;i<nnod;i++) {
		TPZGraphNode *n = &fNodeMap[i];
		if(n) n->DrawSolution(scalind, EMVStyle);
	}
	(fOutFile) << "%RESULT.CASE.STEP.NODAL.VECTOR" << endl;
	(fOutFile) << numvec << endl;
	for(n=0;n<numvec;n++) {
		(fOutFile) << "'" << fVecNames[n] << "' ";
	}
	(fOutFile) << endl;
	(fOutFile) << "%RESULT.CASE.STEP.NODAL.VECTOR.DATA" << endl;
	(fOutFile) << NPoints() << endl;
   nnod = fNodeMap.NElements();
   for(i=0;i<nnod;i++) {
		TPZGraphNode *n = &fNodeMap[i];
		if(n) n->DrawSolution(vecind, EMVStyle);
	}
	if(fNumSteps == fNumCases) (fOutFile) << "%END\n";
}


void TPZMVGraphMesh::SequenceNodes(){
	TPZGraphMesh::SequenceNodes();
   int dim;
   for(dim=0; dim<3; dim++) {
        int nnod = fNodeMap.NElements();
        for(int i=0;i<nnod;i++) {
           TPZGraphNode *n = &fNodeMap[i];
           if(n) n->SetPointNumber(n->FirstPoint()+1);// renumera de 1 para frente
        }                                             // o valor do id do n�
   }
}

void TPZMVGraphMesh::DrawNodes(){

	long nn = 0L;
   int nnod = fNodeMap.NElements();
   int i;
   for(i=0;i<nnod;i++) {
		TPZGraphNode *n = &fNodeMap[i];
		if(n) nn += n->NPoints();
	}
	(fOutFile) << "%NODE" << endl;
	(fOutFile) << nn << endl;

	(fOutFile) << "%NODE.COORD" << endl;
	(fOutFile) << nn << endl;
   for(i=0;i<nnod;i++) {
		TPZGraphNode *n = &fNodeMap[i];
		if(n) n->DrawCo(EMVStyle);
	}
}

void TPZMVGraphMesh::DrawConnectivity(MElementType type) {

	int nel = fElementList.NElements();
	if(!nel) return;
   TPZGraphEl *el = (TPZGraphEl *) fElementList[0];
   int numnodes = el->NConnects();
	(fOutFile) << "%ELEMENT\n";
	int imax = (1<<fResolution);
	(fOutFile) << ((imax*imax)*nel) << endl;
	if(numnodes == 9) {
		(fOutFile) << "%ELEMENT.Q4\n";
   } else {
		(fOutFile) << "%ELEMENT.T3\n";
	}
	(fOutFile) << ((imax*imax)*nel) << endl;
	for(int i=0;i<nel;i++) {
      el = (TPZGraphEl *) fElementList[i];
      if(el) el->Connectivity(EMVStyle);
 	}
}

void TPZMVGraphMesh::DrawSolution(TPZBlock &/*Sol*/) {
    cout << "TPZMVGraphMesh::DrawSolution not Implemented\n";
}

void TPZMVGraphMesh::DrawSolution(char * /*var = 0*/) {
    cout << "TPZMVGraphMesh::DrawSolution not Implemented\n";
}
