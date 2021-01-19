//
//  TPZBuildSBFemHdiv.h
//  PZ
//
//  Created by Karolinne Coelho on 18/01/21.
//
//

#pragma once

#include "TPZBuildSBFem.h"

class TPZBuildSBFemHdiv : public TPZBuildSBFem
{
    //Order of MatIds: fLeftpressure, fRightpressure, fLeftflux, fRightflux;
    TPZManVector<int, 4> fMatIdsHdiv;
    
public:
    
    /// simple constructor
    TPZBuildSBFemHdiv(TPZAutoPointer<TPZGeoMesh> gmesh, int skeletonmatid, std::map<int,int> &matidtranslation) : TPZBuildSBFem(gmesh,skeletonmatid, matidtranslation)
    {
        fMatIdsHdiv.Resize(4);
        fMatIdsHdiv[0] = fSkeletonMatId+1;
        if (fMatIdTranslation[0] == fMatIdsHdiv[0])
        {
            fMatIdsHdiv[0]++;
            std::cout << "TPZBuildSBFemHdiv:: Already using MatId " << fMatIdTranslation[0] << "using MatId" << fMatIdsHdiv[0] << "instead\n";
        }

        for (auto i = 1; i < 4; i++)
        {
            fMatIdsHdiv[i] = fMatIdsHdiv[i-1]+1;
            if (fMatIdTranslation[0] == fMatIdsHdiv[i])
            {
                fMatIdsHdiv[i]++;
                std::cout << "TPZBuildSBFemHdiv:: Already using MatId " << fMatIdTranslation[0] << "using MatId" << fMatIdsHdiv[i] << "instead\n";
            }
        }
    }

    void CreateExternalElements();

    void BuildComputationMeshHdiv(TPZCompMesh &cmesh);

    // void StandardConfigurationHdiv();
    
};