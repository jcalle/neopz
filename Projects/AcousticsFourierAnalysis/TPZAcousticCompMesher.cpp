#include <TPZVTKGeoMesh.h>
#include "TPZAcousticCompMesher.h"
#ifdef STATE_COMPLEX
#include "TPZMatAcousticsFourier.h"
#else
#include "TPZMatAcousticsTransient.h"
#endif
#include "pzbndcond.h"

TPZAcousticCompMesher::TPZAcousticCompMesher(TPZAcousticGeoMesher *geoMesh, bool isAxisymmetric,
        bool filterBoundEqs) :
fGeoMesh(geoMesh) , fCmesh(nullptr), fIsAxisymetric(isAxisymmetric) , fNeqOriginal(-1), fNeqReduced(-1),
fFilterBoundaryEquations(filterBoundEqs){
#ifdef PZDEBUG
    //TODO::Verify if everything is ok
#endif
    fActiveEquations.Resize(0);
}

void
TPZAcousticCompMesher::PrintMesh(const std::string &fileName, const std::string &prefix, bool printVTK, bool printTxt) {

    const std::string meshFileName = prefix + fileName;
    if(printVTK){
        const std::string vtkFile = meshFileName + ".vtk";
        std::ofstream outVTK(vtkFile.c_str());
        TPZVTKGeoMesh::PrintCMeshVTK(fCmesh, outVTK, true);
        outVTK.close();
    }
    if(printTxt){
        const std::string txtFile = meshFileName + ".txt";
        std::ofstream outTXT(txtFile.c_str());
        fCmesh->Print(outTXT);
        outTXT.close();
    }
}



void TPZAcousticCompMesher::CreateFourierMesh(const int &porder) {
    #ifdef STATE_COMPLEX
    fMeshType = ECompMeshTypes::freqDomain;
    TPZAcousticCompMesher::CreateCompMesh<TPZMatAcousticsFourier>(porder);
    #else
    PZError<<"You are trying to perform a frequency-domain analysis without complex numbers (STATE!=complex<double>).";
    PZError<<"Aborting now..."<<std::endl;
    DebugStop();
    #endif
}

void TPZAcousticCompMesher::CreateTransientMesh(const int &porder) {
    #ifdef STATE_REAL
    fMeshType = ECompMeshTypes::timeDomain;
    TPZAcousticCompMesher::CreateCompMesh<TPZMatAcousticsTransient>(porder);
    #else
    PZError<<"You are trying to perform a time-domain analysis with complex numbers (STATE=complex<double>).";
    PZError<<"You are not wrong but we need to implement it. Sorry."<<std::endl;
    PZError<<"Aborting now..."<<std::endl;
    DebugStop();
    #endif
}

void TPZAcousticCompMesher::FilterBoundaryEquations(TPZVec<int64_t> &activeEquations, int64_t &neq,
                                                    int64_t &neqOriginal) {
    neqOriginal = fCmesh->NEquations();
    neq = 0;

    if(fFilterBoundaryEquations){
        std::cout << "Filtering boundary equations..." << std::endl;
        TPZManVector<int64_t, 1000> allConnects;
        std::set<int64_t> boundConnects;

        for (int iel = 0; iel < fCmesh->NElements(); iel++) {
            TPZCompEl *cel = fCmesh->ElementVec()[iel];
            if (cel == NULL) {
                continue;
            }
            if (cel->Reference() == NULL) {
                continue;
            }
            TPZBndCond *mat = dynamic_cast<TPZBndCond *>(fCmesh->MaterialVec()[cel->Reference()->MaterialId()]);
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

        for (int iCon = 0; iCon < fCmesh->NConnects(); iCon++) {
            if (boundConnects.find(iCon) == boundConnects.end()) {
                TPZConnect &con = fCmesh->ConnectVec()[iCon];
                int seqnum = con.SequenceNumber();
                int pos = fCmesh->Block().Position(seqnum);
                int blocksize = fCmesh->Block().Size(seqnum);
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

    return;//TODO: aprender a filtrar eqs
}