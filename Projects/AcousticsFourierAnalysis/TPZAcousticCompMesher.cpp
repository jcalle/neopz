#include <TPZVTKGeoMesh.h>
#include "TPZAcousticCompMesher.h"
#ifdef STATE_COMPLEX
#include "TPZMatAcousticsFourier.h"
#else
#include "TPZMatAcousticsTransient.h"
#endif
#include "pzbndcond.h"

TPZAcousticCompMesher::TPZAcousticCompMesher(TPZAcousticGeoMesher *geoMesh, bool isAxisymmetric) :
fGeoMesh(geoMesh) , fCmesh(nullptr), fIsAxisymetric(isAxisymmetric){
#ifdef PZDEBUG
    //TODO::Verify if everything is ok
#endif
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