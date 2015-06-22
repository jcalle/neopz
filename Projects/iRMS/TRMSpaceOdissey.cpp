//
//  TRMSpaceOdissey.cpp
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#include "TRMSpaceOdissey.h"
#include "TRMFlowConstants.h"

#include "TRMMixedDarcy.h"
#include "TPZMatLaplacian.h"
#include "pzbndcond.h"

#include "tpzgeoelrefpattern.h"
#include "TPZRefPatternTools.h"
#include "tpzhierarquicalgrid.h"
#include "pzgeopoint.h"
#include "TRMSimworxMeshGenerator.h"


static void CreateExampleRawData(TRMRawData &data)
{
    data.fLw = 500.;
    data.fHasLiner = false; //AQUINATHAN esta false para gerar uma malha sem os refinamentos do meio que geram hangnodes
    data.fHasCasing = false; //AQUINATHAN esta false para gerar uma malha sem os refinamentos do meio que geram hangnodes
    
    data.fReservoirWidth = 500.;
    data.fReservoirLength = 1000.;
    data.fReservoirHeight = 50.;
    data.fProdVertPosition = 25;
    data.fWellDiam = 0.2159;
}


/** @brief Default constructor */
TRMSpaceOdissey::TRMSpaceOdissey() : fMeshType(TRMSpaceOdissey::EBox)
{
    
}

/** @brief Default desconstructor */
TRMSpaceOdissey::~TRMSpaceOdissey(){
    
}


/** @brief Create a Hdiv computational mesh Hdiv */
void TRMSpaceOdissey::CreateFluxCmesh(){
    
    if(!fGeoMesh)
    {
        std::cout<< "Geometric mesh doesn't exist" << std::endl;
        DebugStop();
    }
    
    int dim = 3;
    int qorder = 1;
    
    const int typeFlux = 1, typePressure = 0;
    TPZFMatrix<STATE> val1(3,2,0.), val2(3,1,0.);
    
    // Malha computacional
    fFluxCmesh = new TPZCompMesh(fGeoMesh);
    
    TRMMixedDarcy * mat = new TRMMixedDarcy(_ReservMatId);
    fFluxCmesh->InsertMaterialObject(mat);
    
    // Bc N
    TPZBndCond * bcN = mat->CreateBC(mat, _ReservoirInletPressure, typeFlux, val1, val2);
    fFluxCmesh->InsertMaterialObject(bcN);
    
    // Bc S
    TPZBndCond * bcS = mat->CreateBC(mat, _ReservoirOutletPressure, typePressure, val1, val2);
    fFluxCmesh->InsertMaterialObject(bcS);
    
    // Bc E
    TPZBndCond * bcE = mat->CreateBC(mat, _ReservoirNonFluxBoundary, typeFlux, val1, val2);
    fFluxCmesh->InsertMaterialObject(bcE);
    
    // Bc W
    TPZBndCond * bcW = mat->CreateBC(mat, _ReservoirNonFluxBoundary, typeFlux, val1, val2);
    fFluxCmesh->InsertMaterialObject(bcW);
    
    // Bc B
    TPZBndCond * bcB = mat->CreateBC(mat, _ReservoirNonFluxBoundary, typeFlux, val1, val2);
    fFluxCmesh->InsertMaterialObject(bcB);
    
    // Bc T
    TPZBndCond * bcT = mat->CreateBC(mat, _ReservoirNonFluxBoundary, typeFlux, val1, val2);
    fFluxCmesh->InsertMaterialObject(bcT);
    
    // Setando Hdiv
    fFluxCmesh->SetDimModel(dim);
    fFluxCmesh->SetDefaultOrder(qorder);
    fFluxCmesh->SetAllCreateFunctionsHDiv();
    fFluxCmesh->AutoBuild();
    
    
#ifdef DEBUG
    std::ofstream out("CmeshFlux.txt");
    fFluxCmesh->Print(out);
#endif
    
}

/** @brief Create a Discontinuous computational mesh L2 */
void TRMSpaceOdissey::CreatePressureCmesh(){
    if(!fGeoMesh)
    {
        std::cout<< "Geometric mesh doesn't exist" << std::endl;
        DebugStop();
    }
    
    int dim = 3;
    int porder = 1;
    
    const int typeFlux = 1, typePressure = 0;
    TPZFMatrix<STATE> val1(3,2,0.), val2(3,1,0.);
    
    // Malha computacional
    fPressureCmesh = new TPZCompMesh(fGeoMesh);
    
    TRMMixedDarcy * mat = new TRMMixedDarcy(_ReservMatId);
    fPressureCmesh->InsertMaterialObject(mat);
    
    // Bc N
    TPZBndCond * bcN = mat->CreateBC(mat, _ReservoirInletPressure, typeFlux, val1, val2);
    fPressureCmesh->InsertMaterialObject(bcN);
    
    // Bc S
    TPZBndCond * bcS = mat->CreateBC(mat, _ReservoirOutletPressure, typePressure, val1, val2);
    fPressureCmesh->InsertMaterialObject(bcS);
    
    // Bc E
    TPZBndCond * bcE = mat->CreateBC(mat, _ReservoirNonFluxBoundary, typeFlux, val1, val2);
    fPressureCmesh->InsertMaterialObject(bcE);
    
    // Bc W
    TPZBndCond * bcW = mat->CreateBC(mat, _ReservoirNonFluxBoundary, typeFlux, val1, val2);
    fPressureCmesh->InsertMaterialObject(bcW);
    
    // Bc B
    TPZBndCond * bcB = mat->CreateBC(mat, _ReservoirNonFluxBoundary, typeFlux, val1, val2);
    fPressureCmesh->InsertMaterialObject(bcB);
    
    // Bc T
    TPZBndCond * bcT = mat->CreateBC(mat, _ReservoirNonFluxBoundary, typeFlux, val1, val2);
    fPressureCmesh->InsertMaterialObject(bcT);
    
    // Setando L2
    fPressureCmesh->SetDimModel(dim);
    fPressureCmesh->SetDefaultOrder(porder);
    
    fPressureCmesh->SetAllCreateFunctionsContinuous();
    fPressureCmesh->ApproxSpace().CreateDisconnectedElements(false);
    fPressureCmesh->AutoBuild();
    
    fPressureCmesh->AdjustBoundaryElements();
    fPressureCmesh->CleanUpUnconnectedNodes();
    
    int ncon = fPressureCmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = fPressureCmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    
    
#ifdef DEBUG
    std::ofstream out("CmeshPress.txt");
    fPressureCmesh->Print(out);
#endif
    
}

/** @brief Create a Mixed computational mesh Hdiv-L2 */
void TRMSpaceOdissey::CreateMixedCmesh(){
    if(!fGeoMesh)
    {
        std::cout<< "Geometric mesh doesn't exist" << std::endl;
        DebugStop();
    }
    
    int dim = 3;
    
    const int typeFlux = 1, typePressure = 0;
    TPZFMatrix<STATE> val1(1,2,0.), val2(1,1,0.);
    
    // Malha computacional
    fMixedFluxPressureCmesh = new TPZCompMesh(fGeoMesh);
    
    // Material medio poroso
    TRMMixedDarcy * mat = new TRMMixedDarcy(_ReservMatId);
    fMixedFluxPressureCmesh->InsertMaterialObject(mat);
    
    
    // Bc N
    val2(0,0) = 20.0;
    TPZBndCond * bcN = mat->CreateBC(mat, _ReservoirInletPressure, typePressure, val1, val2);
    
    // Bc S
    val2(0,0) = 10.0;
    TPZBndCond * bcS = mat->CreateBC(mat, _ReservoirOutletPressure, typePressure, val1, val2);
    
    // Bc E
    val2(0,0) = 0.0;
    TPZBndCond * bcE = mat->CreateBC(mat, _ReservoirNonFluxBoundary, typeFlux, val1, val2);
    
    // Bc W
    val2(0,0) = 0.0;
    TPZBndCond * bcW = mat->CreateBC(mat, _ReservoirNonFluxBoundary, typeFlux, val1, val2);
    
    // Bc B
    val2(0,0) = 0.0;
    TPZBndCond * bcB = mat->CreateBC(mat, _ReservoirNonFluxBoundary, typeFlux, val1, val2);
    
    // Bc T
    val2(0,0) = 0.0;
    TPZBndCond * bcT = mat->CreateBC(mat, _ReservoirNonFluxBoundary, typeFlux, val1, val2);
    
    
    fMixedFluxPressureCmesh->InsertMaterialObject(bcN);
    fMixedFluxPressureCmesh->InsertMaterialObject(bcS);
    fMixedFluxPressureCmesh->InsertMaterialObject(bcE);
    fMixedFluxPressureCmesh->InsertMaterialObject(bcW);
    fMixedFluxPressureCmesh->InsertMaterialObject(bcB);
    fMixedFluxPressureCmesh->InsertMaterialObject(bcT);
    
    fMixedFluxPressureCmesh->SetDimModel(dim);
    fMixedFluxPressureCmesh->SetAllCreateFunctionsMultiphysicElem();
    fMixedFluxPressureCmesh->AutoBuild();
    
    TPZManVector<TPZCompMesh * ,2> meshvector(2);
    
    this->CreateFluxCmesh();
    this->CreatePressureCmesh();
    
    meshvector[0] = fFluxCmesh.operator->();
    meshvector[1] = fPressureCmesh.operator->();
    
    // Transferindo para a multifisica
    TPZBuildMultiphysicsMesh::AddElements(meshvector, fMixedFluxPressureCmesh.operator->());
    TPZBuildMultiphysicsMesh::AddConnects(meshvector, fMixedFluxPressureCmesh.operator->());
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, fMixedFluxPressureCmesh.operator->());
    
    
}

/** @brief Create a H1 computational mesh */
void TRMSpaceOdissey::CreateH1Cmesh()
{
    if(!fGeoMesh)
    {
        std::cout<< "Geometric mesh doesn't exist" << std::endl;
        DebugStop();
    }
    int porder  = 1;
    fH1Cmesh = new TPZCompMesh(fGeoMesh);
    fH1Cmesh->SetDimModel(3);
    
    TPZMatLaplacian *material = new TPZMatLaplacian(_ReservMatId,3);
    fH1Cmesh->InsertMaterialObject(material);

    TPZFNMatrix<1> val1(1,1,0.),val2(1,1,20);
    TPZBndCond *inflow = new TPZBndCond(material,_ConfinementReservBCbottom,0,val1,val2);
    val2(0,0) = 10.;
    TPZBndCond *outflow = new TPZBndCond(material,_ConfinementReservBCtop,0,val1,val2);
    
//    TPZFNMatrix<1> val1(1,1,0.),val2(1,1,20);
//    TPZBndCond *inflow = new TPZBndCond(material,_ReservoirInletPressure,0,val1,val2);
//    val2(0,0) = 10.;
//    TPZBndCond *outflow = new TPZBndCond(material,_ReservoirOutletPressure,0,val1,val2);
    
    fH1Cmesh->InsertMaterialObject(inflow);
    fH1Cmesh->InsertMaterialObject(outflow);
    fH1Cmesh->SetDefaultOrder(porder);
    
    TPZCreateApproximationSpace space;
    space.SetAllCreateFunctionsContinuous();    
    fH1Cmesh->ApproxSpace() = space;
    
    fH1Cmesh->AutoBuild();
    
#ifdef DEBUG
    std::ofstream out("CmeshPressH1.txt");
    fH1Cmesh->Print(out);
#endif
    
}

void TRMSpaceOdissey::PrintGeometry()
{
    //  Print Geometrical Base Mesh
    std::ofstream argument("GeometicMesh.txt");
    fGeoMesh->Print(argument);
    std::ofstream Dummyfile("GeometricMesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(fGeoMesh,Dummyfile, true);
}

/** @brief Create a reservoir-box geometry */
void TRMSpaceOdissey::CreateGeometricBoxMesh(TPZManVector<int,2> dx, TPZManVector<int,2> dy, TPZManVector<int,2> dz){
    
    REAL t=0.0;
    REAL dt;
    int n;
    
    // Creating a 0D element to be extruded
    TPZGeoMesh * GeoMesh0D = new TPZGeoMesh;
    GeoMesh0D->NodeVec().Resize(1);
    TPZGeoNode Node;
    TPZVec<REAL> coors(3,0.0);
    Node.SetCoord(coors);
    Node.SetNodeId(0);
    GeoMesh0D->NodeVec()[0]=Node;
    
    TPZVec<long> Topology(1,0);
    int elid=0;
    
    new TPZGeoElRefPattern < pzgeom::TPZGeoPoint >(elid,Topology,_ReservMatId,*GeoMesh0D);
    GeoMesh0D->BuildConnectivity();
    GeoMesh0D->SetDimension(0);
    
    TPZHierarquicalGrid CreateGridFrom0D(GeoMesh0D);
    TPZAutoPointer<TPZFunction<STATE> > ParFuncX = new TPZDummyFunction<STATE>(ParametricfunctionX);
    CreateGridFrom0D.SetParametricFunction(ParFuncX);
    CreateGridFrom0D.SetFrontBackMatId(_ReservoirInletPressure,_ReservoirOutletPressure);
    
    dt = dx[0];
    n = dx[1];
    // Computing Mesh extruded along the parametric curve Parametricfunction
    TPZGeoMesh * GeoMesh1D = CreateGridFrom0D.ComputeExtrusion(t, dt, n);
    
    TPZHierarquicalGrid CreateGridFrom1D(GeoMesh1D);
    TPZAutoPointer<TPZFunction<STATE> > ParFuncY = new TPZDummyFunction<STATE>(ParametricfunctionY);
    CreateGridFrom1D.SetParametricFunction(ParFuncY);
    CreateGridFrom1D.SetFrontBackMatId(_ReservoirNonFluxBoundary,_ReservoirNonFluxBoundary);
    CreateGridFrom1D.SetTriangleExtrusion();
    
    dt = dy[0];
    n = dy[1];
    // Computing Mesh extruded along the parametric curve Parametricfunction2
    TPZGeoMesh * GeoMesh2D = CreateGridFrom1D.ComputeExtrusion(t, dt, n);
    
    
    TPZHierarquicalGrid CreateGridFrom2D(GeoMesh2D);
    TPZAutoPointer<TPZFunction<STATE> > ParFuncZ = new TPZDummyFunction<STATE>(ParametricfunctionZ);
    CreateGridFrom2D.SetParametricFunction(ParFuncZ);
    CreateGridFrom2D.SetFrontBackMatId(_ReservoirNonFluxBoundary,_ReservoirNonFluxBoundary);
    CreateGridFrom2D.SetTriangleExtrusion();
    
    dt = dz[0];
    n = dz[1];
    // Computing Mesh extruded along the parametric curve Parametricfunction2
    fGeoMesh = CreateGridFrom2D.ComputeExtrusion(t, dt, n);
    
}
void TRMSpaceOdissey::ParametricfunctionX(const TPZVec<STATE> &par, TPZVec<STATE> &X)
{
    X[0] = par[0];
    X[1] = 0.0;
    X[2] = 25.0*sin(0.1*par[0]);
}

void TRMSpaceOdissey::ParametricfunctionY(const TPZVec<STATE> &par, TPZVec<STATE> &X)
{
    X[0] = 0.0;
    X[1] = par[0];
    X[2] = 0.0;
}

void TRMSpaceOdissey::ParametricfunctionZ(const TPZVec<STATE> &par, TPZVec<STATE> &X)
{
    X[0] = 0.0;
    X[1] = 0.0;
    X[2] = par[0];
}



/** @brief Create the reservoir geometry */
void TRMSpaceOdissey::CreateGeometricReservoirMesh(){
    
    gRefDBase.ReadRefPatternDBase("../RefPatterns.rpt");
    TRMRawData rawdata;
    CreateExampleRawData(rawdata);
    TRMSimworxMeshGenerator meshGen;
    fGeoMesh = meshGen.CreateSimworxGeoMesh(rawdata);
    
}

