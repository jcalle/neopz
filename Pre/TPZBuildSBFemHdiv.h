//
//  TPZBuildSBFemHdiv.h
//  PZ
//
//  Created by Karolinne Coelho on 18/01/21.
//
//

#pragma once

#include "TPZBuildSBFem.h"
#include "TPZMultiphysicsCompMesh.h"

class TPZBuildSBFemHdiv : public TPZBuildSBFem
{
    //Order of MatIds: fLeftpressure, fRightpressure, fLeftflux, fRightflux;
    // TPZManVector<int, 4> fMatIdsHdiv;
    int fLeftfluxMatId, fRightfluxMatId;

    int fLeftpressureMatId, fRightpressureMatId;

    int fInterfaceMatId;
    
public:
    
    /// simple constructor
    TPZBuildSBFemHdiv(TPZAutoPointer<TPZGeoMesh> &gmesh, int skeletonmatid, std::map<int,int> &matidtranslation) : TPZBuildSBFem(gmesh,skeletonmatid, matidtranslation)
    {
        fLeftpressureMatId = skeletonmatid+1;
        fRightpressureMatId = fLeftpressureMatId+1;
        fLeftfluxMatId = fRightpressureMatId+1;
        fRightfluxMatId = fLeftfluxMatId+1;
        fInterfaceMatId = fRightfluxMatId+1;
    }

    void CreateExternalElements(TPZGeoMesh * gmesh, std::set<int> & matids);

    void BuildMultiphysicsCompMesh(TPZCompMesh &cmesh);

    void BuildMultiphysicsCompMeshfromSkeleton(TPZCompMesh &cmesh);

    void CreateSBFemDiscontinuousElements(TPZCompMesh &cmeshpressure);

    void CreateSBFemHdivElements(TPZCompMesh &cmeshf);

    void CreateSBFemMultiphysicsElements(TPZManVector<TPZCompMesh *,2> &cmeshvec, TPZMultiphysicsCompMesh & cmeshm);

    void CreateSBFemInterfaceElementGroups(TPZCompMesh & cmeshm);

    void CreateMultiphysicsElementGroups(TPZMultiphysicsCompMesh * cmesh);

    void AdjustExternalPressureConnectivity();

    void GroupandCondense();
    
};