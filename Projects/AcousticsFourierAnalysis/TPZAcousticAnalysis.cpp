#include "TPZAcousticAnalysis.h"
#include "TPZAcousticCompMesher.h"
#include "TPZSpStructMatrix.h"
#include "pzbndcond.h"



TPZAcousticAnalysis::TPZAcousticAnalysis(TPZAcousticCompMesher *compMesher, const int &nThreads,
        const REAL &deltaT, const int &nTimeSteps,const bool &filter) :
fNThreads(nThreads), fPzAnalysis(compMesher->fCmesh), fCompMesher(compMesher), fNeqReduced(-1),
fFilterBoundaryEquations(filter)
{
    fNTimeSteps = nTimeSteps;
    fDeltaT = deltaT;
    fTotalTime = (fNTimeSteps-1) * fDeltaT;
    TPZSpStructMatrix structMatrix(fCompMesher->fCmesh);
    structMatrix.SetNumThreads(nThreads);
    fNeqOriginal = fCompMesher->fCmesh->NEquations();
    fActiveEquations.Resize(0);
    if(fFilterBoundaryEquations){
        FilterBoundaryEquations(fActiveEquations, fNeqReduced, fNeqOriginal);
        structMatrix.EquationFilter().SetActiveEquations(fActiveEquations);
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

//    const uint solSize = fFilterBoundaryEquations? fCompMesher->fNeqReduced : fCompMesher->fNeqOriginal;
    const uint solSize =  fNeqOriginal;

    TPZFMatrix<STATE> currentSol(solSize,1);
    for(int iTime = 0; iTime < fNTimeSteps; iTime++){
        std::cout<<"\rtime: "<<iTime+1<<" out of "<<fNTimeSteps<<std::flush;
        for(int iPt = 0; iPt < solSize; iPt++){
            currentSol(iPt,0) = fTimeDomainSolution(iPt,iTime);
        }
        fPzAnalysis.LoadSolution(currentSol);
        fPzAnalysis.PostProcess(postProcessResolution);
    }
    std::cout<<std::endl<<" Done!"<<std::endl;
}

void TPZAcousticAnalysis::FilterBoundaryEquations(TPZVec<int64_t> &activeEquations, int64_t &neq,
                                                    int64_t &neqOriginal) {
    TPZCompMesh *cmesh= fCompMesher->fCmesh;
    neqOriginal = cmesh->NEquations();
    neq = 0;

    if(fFilterBoundaryEquations){
        std::cout << "Filtering boundary equations..." << std::endl;
        TPZManVector<int64_t, 1000> allConnects;
        std::set<int64_t> boundConnects;

        for (int iel = 0; iel < cmesh->NElements(); iel++) {
            TPZCompEl *cel = cmesh->ElementVec()[iel];
            if (cel == NULL) {
                continue;
            }
            if (cel->Reference() == NULL) {
                continue;
            }
            TPZBndCond *mat = dynamic_cast<TPZBndCond *>(cmesh->MaterialVec()[cel->Reference()->MaterialId()]);
            if (mat && mat->Type() == 0) {
                std::set<int64_t> boundConnectsEl;
                cel->BuildConnectList(boundConnectsEl);

                for (std::set<int64_t>::iterator iT = boundConnectsEl.begin();
                     iT != boundConnectsEl.end(); iT++) {
                    const int64_t val = *iT;
                    if (boundConnects.find(val) == boundConnects.end()) {
                        boundConnects.insert(val);
                    }
                }
            }
        }

        for (int iCon = 0; iCon < cmesh->NConnects(); iCon++) {
            if (boundConnects.find(iCon) == boundConnects.end()) {
                TPZConnect &con = cmesh->ConnectVec()[iCon];
                int seqnum = con.SequenceNumber();
                int pos = cmesh->Block().Position(seqnum);
                int blocksize = cmesh->Block().Size(seqnum);
                if (blocksize == 0)
                    continue;

                int vs = activeEquations.size();
                activeEquations.Resize(vs + blocksize);
                for (int ieq = 0; ieq < blocksize; ieq++) {
                    activeEquations[vs + ieq] = pos + ieq;
                    neq++;
                }
            }
        }
        std::cout << "# equations(before): " << neqOriginal << std::endl;
        std::cout << "# equations(after): " << neq << std::endl;
    }else{
        neq = neqOriginal;
        std::cout << "Not filtering equations (debug reasons?)" << std::endl;
        std::cout << "# equations: " << neqOriginal << std::endl;
    }

    return;
}