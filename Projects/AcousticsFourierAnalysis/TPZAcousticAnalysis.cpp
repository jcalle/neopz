#include "TPZAcousticAnalysis.h"
#include "TPZAcousticCompMesher.h"
#include "TPZSpStructMatrix.h"

TPZAcousticAnalysis::TPZAcousticAnalysis(TPZAcousticCompMesher *compMesher, const int &nThreads) :
fFilterBoundaryEquations(compMesher->fFilterBoundaryEquations), fNThreads(nThreads),
fPzAnalysis(compMesher->fCmesh), fCompMesher(compMesher)
{
    nTimeSteps = -1;
    TPZSpStructMatrix structMatrix(fCompMesher->fCmesh);
    structMatrix.SetNumThreads(nThreads);
    if(fFilterBoundaryEquations){
        structMatrix.EquationFilter().SetActiveEquations(fCompMesher->fActiveEquations);
    }
    fPzAnalysis.SetStructuralMatrix(structMatrix);
}

void TPZAcousticAnalysis::PostProcess(int vtkResolution, std::string &prefix) {

    TPZCompMesh * cmesh = fCompMesher->fCmesh;

    std::cout<<"Post processing... "<<std::endl;

    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("Pressure");
    std::string plotfile = prefix+"sol";
    plotfile.append(std::to_string(cmesh->NElements()));
    plotfile.append(".vtk");

    fPzAnalysis.DefineGraphMesh(2, scalnames, vecnames,
                       plotfile);  // define malha grafica
    int postProcessResolution = vtkResolution; // define resolucao do pos processamento

    const uint neqOriginal = fCompMesher->fNeqOriginal;
    const uint neq = fCompMesher->fNeqReduced;

    uint solSize = fCompMesher->fFilterBoundaryEquations ? neqOriginal : neq;

    TPZFMatrix<STATE> currentSol(solSize,1);
    for(int iTime = 0; iTime < nTimeSteps; iTime++){
        std::cout<<"\rtime: "<<iTime+1<<" out of "<<nTimeSteps<<std::flush;
        for(int iPt = 0; iPt < solSize; iPt++){
            currentSol(iPt,0) = fTimeDomainSolution(iPt,iTime);
        }
        fPzAnalysis.LoadSolution(currentSol);
        fPzAnalysis.PostProcess(postProcessResolution);
    }
    std::cout<<std::endl<<" Done!"<<std::endl;
}
