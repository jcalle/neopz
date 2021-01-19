//
//  TPZBuildSBFemHdiv.cpp
//  PZ
//
//  Created by Karolinne Coelho on 18/01/21.
//
//

#include "TPZBuildSBFemHdiv.h"
#include "pzgeoelbc.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.tpzbuildsbfem"));
#endif

// standard configuration means each element is a partition and a center node is created
void TPZBuildSBFemHdiv::BuildComputationMeshHdiv(TPZCompMesh &cmesh)
{
    TPZBuildSBFem::BuildComputationMesh(cmesh);

    CreateExternalElements();

    cmesh.SetReference(fGMesh);
}

void TPZBuildSBFemHdiv::CreateExternalElements()
{
    for (auto gel : fGMesh->ElementVec())
    {
        if (!gel || gel->Dimension() != fGMesh->Dimension()-1)
        {
            continue;
        }
        if (fGMesh->Dimension() == 2)
        {
            TPZGeoElBC(gel,2,fMatIdsHdiv[0]); // fLeftpressure
            TPZGeoElBC(gel,2,fMatIdsHdiv[1]); // fRightpressure 
            TPZGeoElBC(gel,2,fMatIdsHdiv[2]); // fLeftflux
            TPZGeoElBC(gel,2,fMatIdsHdiv[3]); // fRightflux
        }
    }
}