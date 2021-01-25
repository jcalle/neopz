//
//  TPZBuildSBFemHdiv.cpp
//  PZ
//
//  Created by Karolinne Coelho on 18/01/21.
//
//

#include "TPZBuildSBFemHdiv.h"

#include "TPZCompElHDivSBFem.h"
#include "TPZNullMaterial.h"
#include "TPZLagrangeMultiplier.h"
#include "pzgeoelbc.h"
#include "pzshapelinear.h"
#include "pzshapetriang.h"
#include "pzshapequad.h"
#include "pzbndcond.h"
#include "TPZVTKGeoMesh.h"

#include "TPZSBFemMultiphysicsElGroup.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.tpzbuildsbfem"));
#endif

using namespace std;

void TPZBuildSBFemHdiv::CreateExternalElements(TPZGeoMesh * gmesh, set<int> & matids)
{
    set<int> matidstarget;
    // Getting the volumetric matids
    for (auto it = fMatIdTranslation.begin(); it!= fMatIdTranslation.end(); it++) {
        int64_t mat = it->second;
        matidstarget.insert(it->second);
    }

    for (auto gel : gmesh->ElementVec())
    {
        if (!gel)
        {
            continue;
        }
        auto it = matidstarget.find(gel->MaterialId());
        if (it == matidstarget.end() )
        {
            continue;
        }

        auto nnodes = gel->NNodes();
        int iside = -1;
        if (nnodes == 4)
        {
            iside = 4;
        }
        else if (nnodes == 6)
        {
            iside = 15;
        }
        else if (nnodes == 8)
        {
            iside = 20;
        }
        else
        {
            DebugStop();
        }

        TPZGeoElSide gelside(gel,iside);
        auto gel1dside = gelside.HasNeighbour(fSkeletonMatId);
        if (!gel1dside)
        {
            // continue;
            TPZGeoElBC(gel, iside, fSkeletonMatId);
        }
        TPZGeoElBC(gel, iside, fLeftpressureMatId);
        TPZGeoElBC(gel, iside, fRightpressureMatId);
        TPZGeoElBC(gel, iside, fLeftfluxMatId);
        TPZGeoElBC(gel, iside, fRightfluxMatId);
    }
}

void TPZBuildSBFemHdiv::BuildMultiphysicsCompMesh(TPZCompMesh &cmesh)
{
    auto cmeshm = dynamic_cast<TPZMultiphysicsCompMesh * >(&cmesh);
    if (!cmeshm)
    {
        DebugStop();
    }
    
    TPZManVector<TPZCompMesh*, 2> cmeshvec = cmeshm->MeshVector();
    auto cmeshflux = cmeshvec[0];
    auto cmeshpressure = cmeshvec[1];

    // create the lower dimensional mesh
    // Creating the external elements - for the hybridization of SBFem-Hdiv element:
    set<int> matids;
    int dim = cmeshflux->Dimension();
    TPZGeoMesh * gmeshflux = cmeshflux->Reference();
    for (auto gel : gmeshflux->ElementVec()) {
        if (!gel) {
            continue;
        }
        if (gel->Dimension() < dim) {
            matids.insert(gel->MaterialId());
        }
    }
    
    cmeshflux->ApproxSpace().SetAllCreateFunctionsHDiv(dim);
    cmeshflux->AutoBuild(matids);

#ifdef PZDEBUG
    ofstream gout0("buildsbfemgmesh0.txt");
    cmeshflux->Reference()->Print(gout0);
#endif

    // cmeshflux->SetReference();

    // Creating geometric elements based on the lower dimensional elements:
    CreateVolumetricElements(*cmeshflux);

    CreateExternalElements(gmeshflux, matids);
    fGMesh = gmeshflux;
#ifdef PZDEBUG
    ofstream outvtk("GeometrySBFEM.vtk");
    TPZVTKGeoMesh vtk;
    vtk.PrintGMeshVTK(fGMesh, outvtk, true);
    ofstream gout("buildsbfemgmesh.txt");
    cmeshflux->Reference()->Print(gout);
#endif

    // Create discontinuous comp elements for the external pressures
    cmeshpressure->SetReference(fGMesh);
    CreateSBFemDiscontinuousElements(*cmeshpressure);

    // Create SBFem Hdiv comp elements
    CreateSBFemHdivElements(*cmeshflux);

    // Updating the multiphysics mesh
    cmeshvec[0] = cmeshflux;
    cmeshvec[1] = cmeshpressure;

    // Updating the geometric mesh:
    // cmeshm->SetReference(fGMesh);
#ifdef PZDEBUG
    ofstream fout("cmeshflux.txt");
    cmeshflux->Print(fout);
    ofstream pout("cmeshpressure.txt");
    cmeshpressure->Print(pout);
#endif

    cmeshflux->SetReference(fGMesh);
    CreateSBFemMultiphysicsElements(cmeshvec, *cmeshm);

    // Creating Element Groups
    CreateSBFemInterfaceElementGroups(*cmeshm);

    // Ajusting the connectivity of these elements
    // AdjustExternalPressureConnectivity();

    // Group and Condense
    // GroupandCondense();
}

void TPZBuildSBFemHdiv::BuildMultiphysicsCompMeshfromSkeleton(TPZCompMesh &cmesh)
{
    DebugStop();
}

void TPZBuildSBFemHdiv::CreateSBFemDiscontinuousElements(TPZCompMesh &cmeshpressure)
{
    auto dim = cmeshpressure.Dimension()-1;
    auto porder = cmeshpressure.GetDefaultOrder();
    auto nstate = cmeshpressure.MaterialVec().begin()->second->NStateVariables();

    cmeshpressure.CleanUp();
    cmeshpressure.SetReference(fGMesh);
    cmeshpressure.SetDefaultOrder(porder);
    cmeshpressure.SetDimModel(dim);
    
    auto matskeleton = new TPZNullMaterial(fSkeletonMatId, dim, nstate);
    cmeshpressure.InsertMaterialObject(matskeleton);

    auto matleft = new TPZNullMaterial(fLeftpressureMatId, dim, nstate);
    cmeshpressure.InsertMaterialObject(matleft);

    auto matright = new TPZNullMaterial(fRightpressureMatId, dim, nstate);
    cmeshpressure.InsertMaterialObject(matright);

    set<int> matids = {fSkeletonMatId, fLeftpressureMatId, fRightpressureMatId};
    cmeshpressure.SetAllCreateFunctionsContinuous();
    cmeshpressure.ApproxSpace().CreateDisconnectedElements(true);
    cmeshpressure.AutoBuild(matids);

    for(auto newnod : cmeshpressure.ConnectVec())
    {
        newnod.SetLagrangeMultiplier(1);
    }
}

void TPZBuildSBFemHdiv::CreateSBFemHdivElements(TPZCompMesh &cmeshflux)
{
    // Getting geometric mesh
    auto gmesh = cmeshflux.Reference();
    auto dim = gmesh->Dimension()-1;

    // Getting the boundary conditions
    set<int> matids, matidstarget;
    for (auto gel : gmesh->ElementVec())
    {
        auto gelmatid = gel->MaterialId();
        if (gelmatid == fSkeletonMatId || gelmatid == fLeftfluxMatId || gelmatid == fRightfluxMatId)
        {
            continue; // won't be set as BCs
        }
        if (gelmatid == fLeftpressureMatId || gelmatid == fRightpressureMatId)
        {
            continue;
        }
        if (gel->Dimension() == dim && gel->MaterialId() != fSkeletonMatId) {
            matids.insert(gel->MaterialId());
        }
    }
    
    // Getting the volumetric matids
    for (auto it = fMatIdTranslation.begin(); it!= fMatIdTranslation.end(); it++) {
        int64_t mat = it->second;
        if (cmeshflux.FindMaterial(mat)) {
            matidstarget.insert(it->second);
        }
    }

    auto nstate = cmeshflux.MaterialVec().begin()->second->NStateVariables();
    auto porder = cmeshflux.GetDefaultOrder();

    // Cleaning all elements of the SBFem-Hdiv
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

    TPZFMatrix<STATE> val1(dim,2,0.), val2(dim,nstate,0.);
    for (auto mId : matids)
    {
        auto bcond = mat->CreateBC(mat, mId, 0, val1, val2);
        cmeshflux.InsertMaterialObject(bcond);
    }
    
    matids.insert(fSkeletonMatId);
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
            iside = 15;
        }
        else if (nnodes == 8)
        {
            iside = 20;
        }
        else
        {
            DebugStop();
        }

        TPZGeoElSide gelside(gel,iside);
        auto gel1dside = gelside.HasNeighbour(fSkeletonMatId);
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

void TPZBuildSBFemHdiv::CreateSBFemMultiphysicsElements(TPZManVector<TPZCompMesh *,2> &cmeshvec, TPZMultiphysicsCompMesh & cmesh)
{
    auto dim = fGMesh->Dimension();
    auto nstate = cmesh.MaterialVec().begin()->second->NStateVariables();

    set<int> matids = {fLeftfluxMatId, fRightpressureMatId, fLeftpressureMatId, fRightpressureMatId};
    for(auto mId : matids)
    {
        auto mat = new TPZNullMaterial(mId,dim,nstate);
        cmesh.InsertMaterialObject(mat);
    }
    {
        auto mat = new TPZLagrangeMultiplier(fInterfaceMatId, dim, nstate);
        cmesh.InsertMaterialObject(mat);
    }

    TPZManVector<int> active(2,1);
    cmesh.BuildMultiphysicsSpace(active, cmeshvec);
    cmesh.LoadReferences();
    cmesh.CleanUpUnconnectedNodes();
}

void TPZBuildSBFemHdiv::CreateSBFemInterfaceElementGroups(TPZCompMesh & cmeshm)
{
    auto gmesh = cmeshm.Reference();
    for (auto gel : gmesh->ElementVec())
    {
        if (!gel || gel->MaterialId() != fLeftfluxMatId)
        {
            continue;
        }
        auto nsides = gel->NSides();
        TPZGeoElSide gelside(gel,nsides-1);
        auto gelsidepr = gelside.HasNeighbour(fLeftpressureMatId);
        if (!gelsidepr)
        {
            DebugStop();
        }

        TPZCompElSide celside = gelside.Reference();
        TPZCompElSide celneigh = gelsidepr.Reference();
        if (!celside || !celneigh) {
            DebugStop();
        }
        TPZGeoElBC gelbc(gelside, fInterfaceMatId);
        int64_t index;
        TPZMultiphysicsInterfaceElement *intf = new TPZMultiphysicsInterfaceElement(cmeshm,gelbc.CreatedElement(),index,celneigh,celside);
    }

    for (auto gel : gmesh->ElementVec())
    {
        if (!gel || gel->MaterialId() != fRightfluxMatId)
        {
            continue;
        }
        auto nsides = gel->NSides();
        TPZGeoElSide gelside(gel,nsides-1);
        auto gelsidepr = gelside.HasNeighbour(fRightpressureMatId);
        if (!gelsidepr)
        {
            DebugStop();
        }

        TPZCompElSide celside = gelside.Reference();
        TPZCompElSide celneigh = gelsidepr.Reference();
        if (!celside || !celneigh) {
            DebugStop();
        }
        TPZGeoElBC gelbc(gelside, fInterfaceMatId);
        int64_t index;
        TPZMultiphysicsInterfaceElement *intf = new TPZMultiphysicsInterfaceElement(cmeshm,gelbc.CreatedElement(),index,celneigh,celside);
    }
}

void TPZBuildSBFemHdiv::AdjustExternalPressureConnectivity()
{
    DebugStop();
}

void TPZBuildSBFemHdiv::GroupandCondense()
{
    DebugStop();
}