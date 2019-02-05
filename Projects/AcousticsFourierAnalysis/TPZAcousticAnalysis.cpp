#include "TPZAcousticAnalysis.h"
#include "TPZAcousticCompMesher.h"
#include "TPZSpStructMatrix.h"

TPZAcousticAnalysis::TPZAcousticAnalysis(TPZAcousticCompMesher *compMesher, const int &nThreads) :
fFilterBoundaryEquations(compMesher->fFilterBoundaryEquations), fNThreads(nThreads),
fPzAnalysis(compMesher->fCmesh), fCompMesher(compMesher)
{
    TPZSpStructMatrix structMatrix(fCompMesher->fCmesh);
    structMatrix.SetNumThreads(nThreads);
    if(fFilterBoundaryEquations){
        structMatrix.EquationFilter().SetActiveEquations(fCompMesher->fActiveEquations);
    }
    fPzAnalysis.SetStructuralMatrix(structMatrix);
}

void TPZAcousticAnalysis::InitializeSolver() {
    fStepSolver.SetDirect(ELU);
    fPzAnalysis.SetSolver(fStepSolver);
}
