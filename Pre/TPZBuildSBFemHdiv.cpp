//
//  TPZBuildSBFemHdiv.cpp
//  PZ
//
//  Created by Karolinne Coelho on 18/01/21.
//
//

#include "TPZBuildSBFemHdiv.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZCompElHDivSBFem.h"
#include "TPZNullMaterial.h"
#include "pzgeoelbc.h"
#include "pzshapelinear.h"
#include "pzshapetriang.h"
#include "pzshapequad.h"
#include "pzbndcond.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.tpzbuildsbfem"));
#endif

using namespace std;

void TPZBuildSBFemHdiv::CreateExternalElements(TPZGeoMesh * gmesh)
{
    for (auto gel : gmesh->ElementVec())
    {
        if (!gel || gel->MaterialId() != fSkeletonMatId)
        {
            continue;
        }
        auto nnodes = gel->NNodes();
        int iside = -1;
        switch (nnodes)
        {
        case 2: // Skeleton is a line
            iside = 2;
            break;
        
        case 3: // Skeleton is composed by triangles
            iside = 6;
            break;

        case 4: // Skeleton is composed by quadrilaterals
            iside = 8;
            break;
        
        default:
            DebugStop();
        }
        TPZGeoElBC(gel, iside, fLeftpressureMatId);
        TPZGeoElBC(gel, iside, fRightpressureMatId);
        TPZGeoElBC(gel, iside, fLeftfluxMatId);
        TPZGeoElBC(gel, iside, fRightfluxMatId);
    }
}

void TPZBuildSBFemHdiv::BuildMultiphysicsCompMesh(TPZCompMesh &cmesh)
{
    TPZMultiphysicsCompMesh * cmeshm = dynamic_cast<TPZMultiphysicsCompMesh * >(&cmesh);
    if (!cmeshm)
    {
        DebugStop();
    }
    
    auto cmeshvec = cmeshm->MeshVector();
    auto cmeshflux = cmeshvec[0];

    // create the lower dimensional mesh
    std::set<int> matids;
    int dim = cmeshflux->Dimension();
    TPZGeoMesh *gmesh = cmeshflux->Reference();
    int64_t nel = gmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (!gel) {
            continue;
        }
        if (gel->Dimension() < dim) {
            matids.insert(gel->MaterialId());
        }
    }
    // create the boundary elements
    cmeshflux->ApproxSpace().SetAllCreateFunctionsHDiv(dim);
    cmeshflux->AutoBuild(matids);

    CreateVolumetricElements(*cmeshvec[0]);
    CreateExternalElements(cmeshvec[0]->Reference());
    
#ifdef PZDEBUG
    ofstream gout("builssbfemgmesh.txt");
    cmeshvec[0]->Reference()->Print(gout);
#endif
    CreateSBFemHdivElements(*cmeshflux);
    // CreateMultiphysicsElementGroups(*cmesh);
}

void TPZBuildSBFemHdiv::BuildMultiphysicsCompMeshfromSkeleton(TPZCompMesh &cmesh)
{
    DebugStop();
}

void TPZBuildSBFemHdiv::CreateSBFemHdivElements(TPZCompMesh &cmeshflux)
{
    // Getting geometric mesh
    auto gmesh = cmeshflux.Reference();
    auto dim = gmesh->Dimension()-1;
    std::set<int> matids, matidstarget;

    // Getting d-1 dimensional matids
    for (auto gel : gmesh->ElementVec())
    {
        if (gel->Dimension() == dim) {
            matids.insert(gel->MaterialId());
        }
    }
    
    // Getting the volumetric matids
    for (std::map<int,int>::iterator it = fMatIdTranslation.begin(); it!= fMatIdTranslation.end(); it++) {
        int64_t mat = it->second;
        if (cmeshflux.FindMaterial(mat)) {
            matidstarget.insert(it->second);
        }
    }

    auto nstate = cmeshflux.MaterialVec()[0]->NStateVariables();
    auto porder = cmeshflux.GetDefaultOrder();
    cmeshflux.CleanUp();

    cmeshflux.SetDefaultOrder(porder);
    cmeshflux.SetDimModel(dim);
    cmeshflux.SetReference(gmesh);

    auto mat = new TPZNullMaterial(fSkeletonMatId, dim, nstate);
    cmeshflux.InsertMaterialObject(mat);
    
    auto matleft = new TPZNullMaterial(fLeftfluxMatId, dim, nstate);
    cmeshflux.InsertMaterialObject(matleft);

    auto matright = new TPZNullMaterial(fRightfluxMatId, dim, nstate);
    cmeshflux.InsertMaterialObject(matright);

    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    for (auto mId : matids)
    {
        auto bcond = mat->CreateBC(mat, mId, 0, val1, val2);
        cmeshflux.InsertMaterialObject(bcond);
    }

    map<int64_t,TPZCompEl *> geltocel;
    for (auto gel : gmesh->ElementVec())
    {
        if (!gel) continue;

        // Finding collapsed element
        auto gelmatid = gel->MaterialId();
        auto it = matidstarget.find(gelmatid);
        if (it == matidstarget.end())
        {
            continue;
        }

        int iside = -1;
        auto nnodes = gel->NNodes();

        if (nnodes == 4)
        {
            iside = 4;
        }
        else if (nnodes == 6)
        {
            DebugStop();
        }
        else if (nnodes == 8)
        {
            DebugStop();
        }
        else
        {
            DebugStop();
        }

        TPZGeoElSide gelside(gel,iside);
        auto gel1dside = gelside.HasNeighbour(matids);
        if (!gelside)
        {
            DebugStop();
        }
        auto gel1d = gel1dside.Element();

        int64_t index;
        auto celhdivc = new TPZCompElHDivSBFem<pzshape::TPZShapeLinear>(cmeshflux, gel1d, gelside, index);
        geltocel[gel1d->Index()] = celhdivc;
    }
    cmeshflux.SetDimModel(dim);
    cmeshflux.ApproxSpace().SetAllCreateFunctionsHDiv(dim);
    cmeshflux.ExpandSolution();
    gmesh->ResetReference();
    
    for (auto gel : gmesh->ElementVec())
    {
        if (!gel) continue;
        if (gel->MaterialId() != fLeftfluxMatId) continue;

        TPZGeoElSide gelside(gel,2);
        auto intfluxside = gelside.HasNeighbour(fSkeletonMatId);
        auto cel  = geltocel[intfluxside.Element()->Index()];

        auto nnodes = gel->NNodes();
        if (nnodes == 2)
        {
            TPZCompElHDivSBFem<pzshape::TPZShapeLinear> * celhdivc = dynamic_cast<TPZCompElHDivSBFem<pzshape::TPZShapeLinear> * >(cel);
            if (!celhdivc)
            {
                DebugStop();
            }
            int64_t index;
            auto hdivboundleft = new TPZCompElHDivBound2<pzshape::TPZShapeLinear>(cmeshflux,gel,index);
            hdivboundleft->SetConnectIndex(0,celhdivc->ConnectIndex(3));

            gel->ResetReference();

            auto rightfluxside = gelside.HasNeighbour(fRightfluxMatId);
            auto gel1dright = rightfluxside.Element();

            auto hdivboundright = new TPZCompElHDivBound2<pzshape::TPZShapeLinear>(cmeshflux,gel1dright,index);
            hdivboundright->SetConnectIndex(0,celhdivc->ConnectIndex(4));
            
            gel1dright->ResetReference();
        }
        else if (nnodes == 3)
        {
            DebugStop();
        }
        else if (nnodes == 4)
        {
            DebugStop();
        }
        else
        {
            DebugStop();
        }
    }
    
    cmeshflux.LoadReferences();
    cmeshflux.Reference()->ResetReference();
    cmeshflux.CleanUpUnconnectedNodes();
    cmeshflux.ExpandSolution();
}

void TPZBuildSBFemHdiv::CreateMultiphysicsElementGroups(TPZCompMesh &cmesh)
{
    auto numgroups = fPartitionCenterNode.size();
    auto groupelementindices(numgroups);
    
    TPZManVector<int64_t> elementgroupindices(numgroups);
    
    for (int64_t el=0; el<numgroups; el++) {
        int64_t index;
        // new TPZSBFemElementGroup(cmesh,index);
        // elementgroupindices[el] = index;
    }
}