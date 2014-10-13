
//  TPZDarcyAnalysis.cpp
//  PZ
//
//  Created by Nathan Shauer and Omar Duran on 9/8/14.
//
//

#include "pzlog.h"
#include "TPZDarcyAnalysis.h"
#include "TPZMatDarcy2dhdiv.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "TPZVTKGeoMesh.h"
#include "TPZReadGIDGrid.h"
#include "pzbndcond.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzanalysis.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzstepsolver.h"
#include "TPZCompElDisc.h"
#include "pzl2projection.h"
#include "TPZMatfrac1dhdiv.h"
#include "pzcompelwithmem.h"
#ifdef USING_BOOST
#include <boost/math/special_functions/erf.hpp>
#endif

#include "TPZFracAnalysis.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.frac"));
#endif

TPZDarcyAnalysis::TPZDarcyAnalysis(TPZAutoPointer<TPZFracData> Data)
{
    fData = Data;
    fmeshvec.Resize(2);
    fmustStop = false;
}


TPZDarcyAnalysis::~TPZDarcyAnalysis()
{
    fData = NULL;
    fmeshvec.Resize(2);
    fgmesh = NULL;
    fcmeshMixed = NULL;
    for (int i = 0; i < 2; i++) {
        fmeshvec[i] = NULL;
    }
    fLastStepRhs.Redim(0, 0);
    fmustStop = false;
    
}

/** @brief Initial pressure field */
void InitialPressure(const TPZVec<REAL> &pt, TPZVec<STATE> &disp)
{
    //    REAL x = pt[0];
    //    REAL y = pt[1];
    disp[0] = 20.0e6;// 20 MPa
}

/** @brief Analytic pressure field */
void PressureAnal(const TPZVec<REAL> &pt, REAL time, TPZVec<STATE> &sol, TPZFMatrix<STATE> &flux)
{
#ifdef USING_BOOST
    REAL x = pt[0], t=time;
    if (time <= 1.0e-8){t=1.0e-8;}
    sol[0]      =   (sqrt((4.0*t)/(M_PI))*exp(-1.0*(x*x)/(4.0*t))) - x*(1.0-boost::math::erf(x/sqrt(4.0*t)));
    flux(0,0)   =   (1.0-boost::math::erf(x/sqrt(4.0*t)));
#endif
}

void TPZDarcyAnalysis::Run()
{
    // Parametros
    const int nel = 1;
    
    
    // Computing vl as the first memory value on each integration point
    TPZFracAnalysis * FracAnalysis = new TPZFracAnalysis(fData);
    REAL vl = FracAnalysis->RunUntilOpen();
    TPZFMatrix<REAL> vlMatrix(1,1,vl);
    
    const REAL lFrac = fData->ElSize();
    fData->SetLfrac(lFrac);
    
    // Solving initial darcy
    
    // Malha geometrica
    fgmesh = CreateGMesh(nel);
    this->PrintGeometricMesh(fgmesh);
    

    //Indcluding the fist ghost element of frac
//    DarcyGmesh(fgmesh);
//    this->PrintGeometricMesh(fgmesh);
//    CreateMultiphysicsMesh(vlMatrix);

    InsertFracsGmesh(fgmesh);
    this->PrintGeometricMesh(fgmesh);
    CreateMultiphysicsMesh(vlMatrix);
    
    // Analysis
    bool mustOptimizeBandwidth = false;
    TPZAnalysis *an = new TPZAnalysis(fcmeshMixed,mustOptimizeBandwidth);
    TPZSkylineNSymStructMatrix skyl(fcmeshMixed);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELU);
    an->SetSolver(step);
    an->SetStructuralMatrix(skyl);
//    SolveSistTransient(an,true);

    
    while (fmustStop == false)
    {
        
        
        bool propagate = SolveSistTransientWithFracture(an);
        if (propagate)
        {
            
            // Novo comprimento de fratura
            REAL newLfrac = fData->Lfrac() + fData->ElSize();
            fData->SetLfrac(newLfrac);
            std::cout << "Lfrac = " << newLfrac << std::endl;
            
            
            // create new elements
            
            // put the acum vl allocated on the new element memory
            
            // 
            
            
        }
    
    }
    
    
    delete an;
    
}


bool TPZDarcyAnalysis::SolveSistTransientWithFracture(TPZAnalysis *an)
{
    
    bool propagate = false;
    int nfracel = this->HowManyFracElement();
    int it = 0;
    while (fmustStop == false && propagate == false) {

        this->SetPressureOnNewElement(an);
        an->Solution().Print("ansol");
        AssembleLastStep(an);
        TPZFMatrix<STATE> lastSol = an->Solution();
    
        if (it == 0)
        { // aqui inicializo chutes iniciais para newton depois da propagacao
        
            if(nfracel == 1){ // esse caso eh para o primeiro elemento da simulacao
                an->Solution().Print("anBefore");
                this->ComputeFirstSolForOneELement(an);
                an->Solution().Print("anAfter");
            }
            else{ // aqui eh quando ha mais de 1 elemento
                
                an->Solution().Print("anBefore");
                this->SetPressureOnLastElement(an);
                an->Solution().Print("anAfter");
                
            }
            
            fData->SetAccumVl(0.); // Zerando accumvl para os proximos
            
        }
        
       
        IterativeProcess(an, std::cout, 50);
        
        an->Solution().Print("Solutionerrado");
//        this->PostProcessVTK(an); //AQUIOMAR APAGAR
        
        const REAL qtip = this->Qtip();
        
        if (qtip < 0.) {
            DebugStop();
            propagate = false;
        }
        else{
            propagate = false;//VerifyIfPropagate(qtip);
        }
        
        
        if (propagate) {
            fData->SetLastQtip(qtip);
            an->Solution() = lastSol;
            an->LoadSolution();
            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fcmeshMixed);
            std::cout << "\n******************************** FRACTURE PROPAGATED ********************************" << std::endl;
        }
        else{
            fData->SetNextTime();
            if (qtip > 0.) {
               // AcceptSolution(an); // updates leak off
            }
            this->PostProcessVTK(an);
        }
        
        REAL peteleco = 1.e-8;
        if( fData->Time() > (fData->TotalTime() - peteleco) )
        {
            fmustStop = true;
        }
        it++;
    }
    return propagate;
}


TPZGeoMesh * TPZDarcyAnalysis::CreateGMesh(const int nel)
{
    std::string dirname = PZSOURCEDIR;
    std::string FileName = dirname;
    
    //  Reading mesh
    std::string GridFileName;
    GridFileName = dirname + "/Projects/Frac1DHdiv/";
    GridFileName += "OilWaterSystemUnit.dump";
//    GridFileName += "BaseGeometryDakeThin.dump";//"FiveSpot.dump";
//    GridFileName += "FiveSpot.dump";//"FiveSpot.dump";
    REAL angle = 0.0*M_PI/4.0;
    
    TPZReadGIDGrid GeometryInfo;
    GeometryInfo.SetfDimensionlessL(1.0);
    TPZGeoMesh * gmesh = GeometryInfo.GeometricGIDMesh(GridFileName);
    
//    // Inserting
//    TPZGeoEl *  Quad = gmesh->ElementVec()[0];
//    TPZGeoElSide gelside(Quad,0);
//    TPZGeoEl *  InletPoint = Quad->CreateBCGeoEl(0,7);
//    long index;
//    TPZVec<long> TopologyPoint(1,InletPoint->Node(0).Id());
//    gmesh->CreateGeoElement(EPoint, TopologyPoint, 7, index);
//    RotateGeomesh(gmesh, angle);
    
    UniformRefinement(gmesh, nel);

    //  Print Geometrical Base Mesh
    this->PrintGeometricMesh(gmesh);
    
    return gmesh;
}

TPZCompMesh * TPZDarcyAnalysis::CreateCMeshFluxHdiv()
{
    const int matId2d = 1, matId1d = 6, bcBottomId = 2, bcRightId = 3, bcTopId = 4, bcLeftId = 5, bcBottomIdAux = 10;
    const int MatINterface = 20;
    const int bcinlet = 7, bcoutlet = 8;
    const int typeFlux = 0, typePressure = 1;
    TPZFMatrix<STATE> val1(3,2,0.), val2(3,1,0.);
    
    // Malha computacional
    TPZCompMesh *cmesh = new TPZCompMesh(fgmesh);
    
    TPZMatDarcy2dhdiv *mat = new TPZMatDarcy2dhdiv(matId2d);
    cmesh->InsertMaterialObject(mat);
    
    // Bc Bottom
    TPZBndCond * bcBottom = mat->CreateBC(mat, bcBottomId, typeFlux, val1, val2);
 	cmesh->InsertMaterialObject(bcBottom);
    
    // Bc Right
    TPZBndCond * bcRight = mat->CreateBC(mat, bcRightId, typePressure, val1, val2);
 	cmesh->InsertMaterialObject(bcRight);
    
    // Bc Top
    TPZBndCond * bcTop = mat->CreateBC(mat, bcTopId, typeFlux, val1, val2);
 	cmesh->InsertMaterialObject(bcTop);
    
    // Bc Left
    TPZBndCond * bcLeft = mat->CreateBC(mat, bcLeftId, typeFlux, val1, val2);
 	cmesh->InsertMaterialObject(bcLeft);
    
    // Bc Left
    TPZBndCond * AuxbcLeft = mat->CreateBC(mat, bcBottomIdAux, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(AuxbcLeft);

    // Interface material
    TPZMatfrac1dhdiv *mat1d = new TPZMatfrac1dhdiv(matId1d);
    cmesh->InsertMaterialObject(mat1d);

    // Material da fratura
    TPZMatfrac1dhdiv *matInterface1d = new TPZMatfrac1dhdiv(MatINterface);
    cmesh->InsertMaterialObject(matInterface1d);
    
    // Condicao de contorno na esquerda
    val2(0,0) = fData->Q();
    TPZBndCond * bcin = mat1d->CreateBC(mat1d, bcinlet, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcin);
    
    // Condicao de contorno na direita
    val2(0,0) = fData->SigmaConf();
    TPZBndCond * bcout = mat1d->CreateBC(mat1d, bcoutlet, typePressure, val1, val2);
    cmesh->InsertMaterialObject(bcout);
    
    
    // Setando Hdiv
    cmesh->SetDimModel(2);
    cmesh->SetDefaultOrder(fData->PorderDarcyFlow());
    cmesh->SetAllCreateFunctionsHDiv();
    
    std::set<int> BuildGroup2D,BuildGroup1D;
    BuildGroup2D.insert(matId2d);
    BuildGroup2D.insert(bcLeftId);
    BuildGroup2D.insert(bcRightId);
    BuildGroup2D.insert(bcTopId);
    BuildGroup2D.insert(bcBottomId);
    BuildGroup2D.insert(bcBottomIdAux);
    cmesh->AutoBuild(BuildGroup2D);

    // Setando H1
    fgmesh->ResetReference();
    cmesh->SetDimModel(1);
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->SetDefaultOrder(fData->PorderFlow());
    BuildGroup1D.insert(matId1d);
    BuildGroup1D.insert(bcinlet);
    BuildGroup1D.insert(bcoutlet);
    cmesh->AutoBuild(BuildGroup1D);
 
    cmesh->SetDimModel(2);
    
    return cmesh;
}

TPZCompMesh * TPZDarcyAnalysis::CreateCMeshPressureL2()
{

    const int matId2d = 1, matId1d = 6, bcBottomId = 2, bcRightId = 3, bcTopId = 4, bcLeftId = 5, bcBottomIdAux = 10;
    const int MatINterface = 20;
    const int bcinlet = 7, bcoutlet = 8;
    const int typeFlux = 0, typePressure = 1;
    TPZFMatrix<STATE> val1(3,2,0.), val2(3,1,0.);
    
    // Malha computacional
    TPZCompMesh *cmesh = new TPZCompMesh(fgmesh);
    
    // Material da fratura
    TPZMatDarcy2dhdiv *mat = new TPZMatDarcy2dhdiv(matId2d);
    cmesh->InsertMaterialObject(mat);
    
    // Bc Bottom
    TPZBndCond * bcBottom = mat->CreateBC(mat, bcBottomId, typeFlux, val1, val2);
 	cmesh->InsertMaterialObject(bcBottom);
    
    // Bc Right
    TPZBndCond * bcRight = mat->CreateBC(mat, bcRightId, typePressure, val1, val2);
 	cmesh->InsertMaterialObject(bcRight);
    
    // Bc Top
    TPZBndCond * bcTop = mat->CreateBC(mat, bcTopId, typeFlux, val1, val2);
 	cmesh->InsertMaterialObject(bcTop);
    
    // Bc Left
    TPZBndCond * bcLeft = mat->CreateBC(mat, bcLeftId, typeFlux, val1, val2);
 	cmesh->InsertMaterialObject(bcLeft);
    
    // Bc Left
    TPZBndCond * AuxbcLeft = mat->CreateBC(mat, bcBottomIdAux, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(AuxbcLeft);
    
    // Material da fratura
    TPZMatfrac1dhdiv *matInterface1d = new TPZMatfrac1dhdiv(MatINterface);
    cmesh->InsertMaterialObject(matInterface1d);
    
    // Material da fratura
    TPZMatfrac1dhdiv *mat1d = new TPZMatfrac1dhdiv(matId1d);
    cmesh->InsertMaterialObject(mat1d);
    
    // Condicao de contorno na esquerda
    val2(0,0) = fData->Q();
    TPZBndCond * bcin = mat1d->CreateBC(mat1d, bcinlet, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcin);
    
    // Condicao de contorno na direita
    val2(0,0) = fData->SigmaConf();
    TPZBndCond * bcout = mat1d->CreateBC(mat1d, bcoutlet, typePressure, val1, val2);
    cmesh->InsertMaterialObject(bcout);
    
    
    // Setando L2
    cmesh->SetDimModel(2);
    cmesh->SetDefaultOrder(fData->PorderDarcyPressure());
    cmesh->SetAllCreateFunctionsDiscontinuous();

    std::set<int> BuildGroup2D,BuildGroup1D;
    BuildGroup2D.insert(matId2d);
    BuildGroup2D.insert(bcLeftId);
    BuildGroup2D.insert(bcRightId);
    BuildGroup2D.insert(bcTopId);
    BuildGroup2D.insert(bcBottomId);
    BuildGroup2D.insert(bcBottomIdAux);
    cmesh->AutoBuild(BuildGroup2D);
    
    fgmesh->ResetReference();
    cmesh->SetDimModel(1);
    cmesh->SetDefaultOrder(fData->PorderPressure());
    BuildGroup1D.insert(matId1d);
    BuildGroup1D.insert(bcinlet);
    BuildGroup1D.insert(bcoutlet);
    cmesh->AutoBuild(BuildGroup1D);
    
    cmesh->SetDimModel(2);
    
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    
    return cmesh;
    
}

TPZCompMesh * TPZDarcyAnalysis::CreateCMeshMixed(TPZFMatrix<REAL> Vl)
{
    // Definicao de ids e tipos
    const int matId2d = 1, matId1d = 6, bcBottomId = 2, bcRightId = 3, bcTopId = 4, bcLeftId = 5, bcBottomIdAux = 10;
    const int MatINterface = 20;    
    const int bcinlet = 7, bcoutlet = 8;
    const int typeFlux = 0, typePressure = 1;
    TPZFMatrix<STATE> val1(3,2,0.), val2(3,1,0.);
    
    // Malha computacional
    TPZCompMesh *cmesh = new TPZCompMesh(fgmesh);
    
    // Material medio poroso
    TPZMatDarcy2dhdiv *mat = new TPZMatDarcy2dhdiv(matId2d);
    mat->SetSimulationData(fData);
    cmesh->InsertMaterialObject(mat);
    
    TPZAutoPointer<TPZFunction<STATE> > TimeDepFExact = new TPZDummyFunction<STATE>(PressureAnal);
    mat->SetTimeDependentFunctionExact(TimeDepFExact);
    
    // Bc Bottom
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    val2(2,0) = 0.0;
    TPZBndCond * bcBottom = mat->CreateBC(mat, bcBottomId, typeFlux, val1, val2);
 	cmesh->InsertMaterialObject(bcBottom);
    
    // Bc Right
    val2(0,0) = 0.0;
    val2(1,0) = 0.0000;
    val2(2,0) = 40.0e6;
    TPZBndCond * bcRight = mat->CreateBC(mat, bcRightId, typePressure, val1, val2);
 	cmesh->InsertMaterialObject(bcRight);
    
    // Bc Top
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    val2(2,0) = 0.0;
    TPZBndCond * bcTop = mat->CreateBC(mat, bcTopId, typeFlux, val1, val2);
 	cmesh->InsertMaterialObject(bcTop);
    
    // Bc Left
    val2(0,0) = 0.0;// Massic flux 5.0 kg/s over 100000 m2
    val2(1,0) = 0.0;
    val2(2,0) = 0.0;
    TPZBndCond * bcLeft = mat->CreateBC(mat, bcLeftId, typeFlux, val1, val2);
 	cmesh->InsertMaterialObject(bcLeft);
    
    // Bc Left
    val2(0,0) = 0.0;// Massic flux 5.0 kg/s over 100000 m2
    val2(1,0) = Vl(0,0)/fData->Time();
    val2(2,0) = 0.0;
    TPZBndCond * AuxbcLeft = mat->CreateBC(mat, bcBottomIdAux, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(AuxbcLeft);
    
    // Material da fratura
    TPZMatfrac1dhdiv *matInterface1d = new TPZMatfrac1dhdiv(MatINterface);
    matInterface1d->SetSimulationData(fData);
    matInterface1d->SetDefaultMem(Vl);
    cmesh->InsertMaterialObject(matInterface1d);
    
    // Material da fratura
    TPZMatfrac1dhdiv *mat1d = new TPZMatfrac1dhdiv(matId1d);
    mat1d->SetSimulationData(fData);
    cmesh->InsertMaterialObject(mat1d);
    mat1d->SetDefaultMem(Vl);
    
    // Condicao de contorno na esquerda
    val2(0,0) = fData->Q();
    TPZBndCond * bcin = mat1d->CreateBC(mat1d, bcinlet, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcin);
    
    // Condicao de contorno na direita
    val2(0,0) = fData->SigmaConf();
    TPZBndCond * bcout = mat1d->CreateBC(mat1d, bcoutlet, typePressure, val1, val2);
    cmesh->InsertMaterialObject(bcout);
    
    // Setando Multifisico
    cmesh->SetDimModel(2);
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    
    std::set<int> BuildGroup2D,BuildGroup1D;
    BuildGroup2D.insert(matId2d);
    BuildGroup2D.insert(bcLeftId);
    BuildGroup2D.insert(bcRightId);
    BuildGroup2D.insert(bcTopId);
    BuildGroup2D.insert(bcBottomId);
    BuildGroup2D.insert(bcBottomIdAux);
    cmesh->AutoBuild(BuildGroup2D);
    
    fgmesh->ResetReference();
    cmesh->SetDimModel(1);
    cmesh->SetAllCreateFunctionsMultiphysicElemWithMem();
    BuildGroup1D.insert(matId1d);
    BuildGroup1D.insert(bcinlet);
    BuildGroup1D.insert(bcoutlet);
    cmesh->AutoBuild(BuildGroup1D);
    
    cmesh->SetDimModel(2);
    
    return cmesh;
}

TPZCompMesh * TPZDarcyAnalysis::L2ProjectionP(TPZGeoMesh *gmesh, int pOrder, TPZVec<STATE> &solini)
{
    /// criar materiais
    int dim = 2;
    TPZL2Projection *material;
    material = new TPZL2Projection(1, dim, 1, solini, pOrder);
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
    TPZAutoPointer<TPZFunction<STATE> > forcef = new TPZDummyFunction<STATE>(InitialPressure);
    material->SetForcingFunction(forcef);
    cmesh->SetAllCreateFunctionsDiscontinuous();
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    cmesh->AutoBuild();
    
//    ///set order total da shape HERE when Triangles are used
//    int nel = cmesh->NElements();
//    for(int i=0; i<nel; i++){
//        TPZCompEl *cel = cmesh->ElementVec()[i];
//        TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc *>(cel);
//        if(celdisc && celdisc->Reference()->Dimension() == cmesh->Dimension())
//        {
//            if(ftriang==true || celdisc->Reference()->Type()==ETriangle) celdisc->SetTotalOrderShape();
//            else celdisc->SetTensorialShape();
//        }
//    }
    
    return cmesh;
    
}


void TPZDarcyAnalysis::CreateInterfaces(TPZCompMesh *cmesh)
{
    fgmesh->ResetReference();
    cmesh->LoadReferences();
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    // Creation of interface elements
    int nel = cmesh->ElementVec().NElements();
    for(int el = 0; el < nel; el++)
    {
        TPZCompEl * compEl = cmesh->ElementVec()[el];
        if(!compEl) continue;
        int index = compEl ->Index();
        if(compEl->Dimension() == cmesh->Dimension())
        {
            TPZMultiphysicsElement * InterpEl = dynamic_cast<TPZMultiphysicsElement *>(cmesh->ElementVec()[index]);
            if(!InterpEl) continue;
            InterpEl->CreateInterfaces();
        }

    }
    
    cmesh->SetAllCreateFunctionsMultiphysicElemWithMem();
    cmesh->SetDefaultOrder(fData->PorderPressure());
    // Creation of interface elements for frac coupling
    for(int el = 0; el < nel; el++)
    {
        TPZCompEl * compEl = cmesh->ElementVec()[el];
        if(!compEl) continue;
        int index = compEl ->Index();
        if(compEl->Reference()->MaterialId() == 6)
        {
            TPZGeoEl * gel = compEl->Reference();
            int side = gel->NSides()-1;
            
            TPZStack < TPZCompElSide > neigh;
            TPZGeoEl * GeoLinearFrac = gel->CreateBCGeoEl(side, 20);
            TPZGeoElSide gelside(gel,side);
            TPZCompElSide Linearside(compEl,side);
            
            gelside.EqualLevelCompElementList(neigh, 0, 0);
            
            if (neigh.size() == 0) {
                continue;
            }
            long gelindex;
            new TPZCompElWithMem <TPZMultiphysicsInterfaceElement >(*cmesh, GeoLinearFrac, gelindex, neigh[0], Linearside);
            
        }
    }
    
}

void TPZDarcyAnalysis::IterativeProcess(TPZAnalysis *an, std::ostream &out, int numiter)
{
	int iter = 0;
	REAL error = 1.e10;
    const REAL tol = 1.e-6; // because the SI unit system
    
    fData->SetCurrentState();
	int numeq = an->Mesh()->NEquations();
	
	TPZFMatrix<STATE> Uatk0(an->Solution());
    TPZFMatrix<STATE> Uatk(Uatk0),DeltaU(Uatk0);
	if(Uatk0.Rows() != numeq) Uatk0.Redim(numeq,1);
    
    an->Assemble();
    an->Rhs() += fLastStepRhs;
    an->Rhs() *= -1.0; //- [R(U0)];
    
    TPZAutoPointer< TPZMatrix<REAL> > matK; // getting X(Uatn)

   
    
//    TPZSkylineNSymStructMatrix skyl(fcmeshMixed);
//    std::set<int> matids; // to be computed
//    matids.insert(6);
//    matids.insert(7);
//    matids.insert(8);
//    matids.insert(20);
//    skyl.SetMaterialIds(matids);
//    TPZFMatrix<STATE> rhsfrac;
//    TPZAutoPointer<TPZGuiInterface> gui = new TPZGuiInterface;
//    TPZAutoPointer<TPZMatrix<STATE> > matfrac = skyl.CreateAssemble(rhsfrac, gui);
//    
//    
//#ifdef LOG4CXX
//    if(logger->isDebugEnabled())
//    {
//        std::stringstream sout;
//        matfrac->Print("matfrac = ", sout,EMathematicaInput);
//        rhsfrac.Print("rhsfrac = ", sout,EMathematicaInput);
//        LOGPZ_DEBUG(logger,sout.str())
//    }
//#endif
    
	while(error > tol && iter < numiter) {
		
#ifdef LOG4CXX
        if(logger->isDebugEnabled())
        {
            std::stringstream sout;
            matK=an->Solver().Matrix();
            matK->Print("matK = ", sout,EMathematicaInput);
            LOGPZ_DEBUG(logger,sout.str())
        }
#endif
        
        // Computing Uatk = Uatn + DeltaU;
		an->Solve();
        DeltaU= an->Solution();
        Uatk = Uatk0 + DeltaU;
        
        //Computing ||DeltaU||
		REAL NormOfDeltaU = Norm(DeltaU);
        
#ifdef LOG4CXX
        if(logger->isDebugEnabled())
        {
            std::stringstream sout;
            DeltaU.Print("DeltaU = ", sout,EMathematicaInput);
            Uatk.Print("Uatk = ", sout,EMathematicaInput);
            LOGPZ_DEBUG(logger,sout.str())
        }
#endif

        
        an->LoadSolution(Uatk); // Loading Uatk
        TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fcmeshMixed);
        an->Assemble();

#ifdef LOG4CXX
        if(logger->isDebugEnabled())
        {
            std::stringstream sout;
            fLastStepRhs.Print("ResAtn = ", sout,EMathematicaInput);
            an->Rhs().Print("Respone = ", sout,EMathematicaInput);            
            LOGPZ_DEBUG(logger,sout.str())
        }
#endif
        
        an->Rhs() += fLastStepRhs;
        
#ifdef LOG4CXX
        if(logger->isDebugEnabled())
        {
            std::stringstream sout;
            an->Rhs().Print("Res = ", sout,EMathematicaInput);
            LOGPZ_DEBUG(logger,sout.str())
        }
#endif
        // Computing ||[R(Uatk)]||
        double ResidualNorm = Norm(an->Rhs());
		double norm = ResidualNorm;
		out << "Iteration n : " << (iter+1) << " : norms ||DeltaU|| e ||[R(Uatk)]|| : " << NormOfDeltaU << " / " << ResidualNorm << std::endl;
        
		if(norm < tol /*|| NormResLambda < tol*/) {
			out << "\nTolerance at n : " << (iter+1) << std::endl;
			out << "\n\nNorm ||DeltaU||  : " << NormOfDeltaU << std::endl;
            out << "\n\nNorm ||[R(Uatk)]||  : " << ResidualNorm << std::endl;
            
		}
        else if( (norm - error) > 1.e-9 ) {
            out << "\nMethod with wrong direction\n";
        }
        
		error = norm;
		iter++;
        Uatk0 = Uatk;
		out.flush();
	}
    
    if (error > tol) {
        DebugStop(); // Something is very wrong
    }
    
}

void TPZDarcyAnalysis::PrintGeometricMesh(TPZGeoMesh * gmesh)
{
    std::ofstream argument("GeometicMesh.txt");
    gmesh->Print(argument);
    std::ofstream Dummyfile("GeometricMesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh,Dummyfile, true);
}


void TPZDarcyAnalysis::PostProcessVTK(TPZAnalysis *an)
{
    const int dim = 2;
    int div =2;
    TPZStack<std::string> scalnames, vecnames;
    std::string plotfile = "2DMixedDarcy.vtk";
    scalnames.Push("Pressure");
    scalnames.Push("PressureExact");
    vecnames.Push("MassVelocity");
    vecnames.Push("MassVelocityExact");
    an->DefineGraphMesh(dim, scalnames, vecnames, fData->PostProcessFileName());
    an->PostProcess(div,dim);
}


void TPZDarcyAnalysis::AssembleLastStep(TPZAnalysis *an)
{
    fData->SetLastState();
    an->Assemble();
    fLastStepRhs = an->Rhs();
    an->Rhs().Print("rhs");
}

void TPZDarcyAnalysis::SolveSyst(TPZAnalysis &an, TPZCompMesh *Cmesh)
{
    
    TPZSkylineStructMatrix full(Cmesh);
    an.SetStructuralMatrix(full);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    //step.SetDirect(ELU);
    an.SetSolver(step);
    an.Run();

}

bool TPZDarcyAnalysis::SolveSistTransient(TPZAnalysis *an, bool initial)
{

    this->PostProcessVTK(an);

    while (fmustStop == false) {
        
        AssembleLastStep(an);
        IterativeProcess(an, std::cout);
        fData->SetNextTime();
        this->PostProcessVTK(an);
        
        if (initial) {
            fmustStop=true;
        }
        
        REAL tricky = 1.E-8;
        if( fData->Time() > (fData->TotalTime() - tricky) )
        {
            fmustStop = true;
        }
    }
}



void TPZDarcyAnalysis::UniformRefinement(TPZGeoMesh *gMesh, int nh)
{
    for ( int ref = 0; ref < nh; ref++ ){
        TPZVec<TPZGeoEl *> filhos;
        long n = gMesh->NElements();
        for ( long i = 0; i < n; i++ ){
            TPZGeoEl * gel = gMesh->ElementVec() [i];
            if (gel->Dimension() == 2 || gel->Dimension() == 1) gel->Divide (filhos);
        }//for i
    }//ref
}

void TPZDarcyAnalysis::RotateGeomesh(TPZGeoMesh *gmesh, REAL CounterClockwiseAngle)
{
    REAL theta = CounterClockwiseAngle;
    // It represents a 3D rotation around the z axis.
    TPZFMatrix<STATE> RotationMatrix(3,3,0.0);
    RotationMatrix(0,0) =   +cos(theta);
    RotationMatrix(0,1) =   -sin(theta);
    RotationMatrix(1,0) =   +sin(theta);
    RotationMatrix(1,1) =   +cos(theta);
    RotationMatrix(2,2) = 1.0;
    TPZVec<STATE> iCoords(3,0.0);
    TPZVec<STATE> iCoordsRotated(3,0.0);
    
    RotationMatrix.Print("Rotation = ");
    
    int NumberofGeoNodes = gmesh->NNodes();
    for (int inode = 0; inode < NumberofGeoNodes; inode++)
    {
        TPZGeoNode GeoNode = gmesh->NodeVec()[inode];
        GeoNode.GetCoordinates(iCoords);
        // Apply rotation
        iCoordsRotated[0] = RotationMatrix(0,0)*iCoords[0]+RotationMatrix(0,1)*iCoords[1]+RotationMatrix(0,2)*iCoords[2];
        iCoordsRotated[1] = RotationMatrix(1,0)*iCoords[0]+RotationMatrix(1,1)*iCoords[1]+RotationMatrix(1,2)*iCoords[2];
        iCoordsRotated[2] = RotationMatrix(2,0)*iCoords[0]+RotationMatrix(2,1)*iCoords[1]+RotationMatrix(2,2)*iCoords[2];
        GeoNode.SetCoord(iCoordsRotated);
        gmesh->NodeVec()[inode] = GeoNode;
    }
    
}

void TPZDarcyAnalysis::InsertFracsGmesh(TPZGeoMesh * gmesh)
{
    int nelements = fgmesh->NElements();
    for (int iel = 0 ; iel < nelements; iel++)
    {
        TPZGeoEl *igel = fgmesh->ElementVec()[iel];
        if (igel->HasSubElement()) {
            continue;
        }
        if(igel->MaterialId()== 2)
        {
            TPZGeoElSide igelside(igel,0);
            TPZGeoElSide neigh = igelside.Neighbour();
            
            while (igelside != neigh) {
                if (neigh.Element()->Dimension() == 0 || (neigh.Element()->MaterialId()==5 && neigh.Element()->Dimension() == 1)) {
                    break;
                }
                neigh = neigh.Neighbour();
            }
            if (igelside == neigh) {
                continue;
            }
            
            TPZGeoEl *Neigel = neigh.Element();
            if (Neigel->MaterialId()==5) {
                igel->SetMaterialId(6);
                long rightnode=igel->NodeIndex(1);
                TPZVec<long> topopoint(1,rightnode);
                long OutLetindex;
                fgmesh->CreateGeoElement(EPoint, topopoint, 8, OutLetindex);
                
                long leftnode=igel->NodeIndex(0);
                topopoint=leftnode;
                fgmesh->CreateGeoElement(EPoint, topopoint, 7, OutLetindex);
                
                fgmesh->ResetReference();
                fgmesh->BuildConnectivity();
                break;
            }else if (Neigel->MaterialId()==8)
            {
                igel->SetMaterialId(6);
                long rightnode=igel->NodeIndex(1);
                Neigel->SetNodeIndex(0,rightnode);
                fgmesh->ResetReference();
                fgmesh->BuildConnectivity();
                break;
            }
            
        }
    }
    
    //---------  Elemento de interface  --------------------------------------------
    gmesh->AddInterfaceMaterial(1,6,20);
    gmesh->AddInterfaceMaterial(6,1,20);
    gmesh->BuildConnectivity();
    
}

void TPZDarcyAnalysis::DarcyGmesh(TPZGeoMesh * gmesh)
{
    int nelements = fgmesh->NElements();
    for (int iel = 0 ; iel < nelements; iel++)
    {
        TPZGeoEl *igel = fgmesh->ElementVec()[iel];
        if (igel->HasSubElement()) {
            continue;
        }
        if(igel->MaterialId()== 2)
        {
            TPZGeoElSide igelside(igel,0);
            TPZGeoElSide neigh = igelside.Neighbour();
            
            while (igelside != neigh) {
                if (neigh.Element()->MaterialId()==5 && neigh.Element()->Dimension() == 1) {
                    break;
                }
                neigh = neigh.Neighbour();
            }
            if (igelside == neigh) {
                continue;
            }
            igel->SetMaterialId(10);
        }
    }
    
}

// Computing the flux on the right tip of the fracture

REAL TPZDarcyAnalysis::Qtip()
{
    fgmesh->ResetReference();
    fcmeshMixed->LoadReferences();
    const int bcOfTip = 6;
    const long nel = fcmeshMixed->NElements();
    TPZCompEl *cel = NULL;
    for (long iel = 0; iel < nel; iel++) {
        cel = fcmeshMixed->Element(iel);
        if (!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if (gel->Dimension() != 1) {
            continue;
        }
        TPZGeoElSide gelside(gel,1);
        TPZGeoElSide neigh = gelside.Neighbour();
        while (neigh != gelside) {
            if (neigh.Element()->MaterialId() == bcOfTip) {
                break;
            }
            neigh = neigh.Neighbour();
        }
        if (neigh == gelside) {
            continue;
        }
        break;
    }
    
    
    
    //Here mat 6 means fracture elements
    TPZMaterial *Material = fcmeshMixed->FindMaterial(bcOfTip);
    TPZMatfrac1dhdiv * MaterialOfFract = dynamic_cast<TPZMatfrac1dhdiv *>(Material);

    
    TPZVec<REAL> qsi(3,1.), sol(1,0.);
    const int varQ = MaterialOfFract->VariableIndex("Flow");
    cel->Solution(qsi, varQ, sol);
    const REAL qTip = sol[0];
    std::cout << "\nqtip = " << qTip << std::endl;
    
    return qTip;
}


void TPZDarcyAnalysis::SetPressureOnLastElement(TPZAnalysis *an)
{
    fgmesh->ResetReference();
    fmeshvec[1]->LoadReferences();
    
    long nfracel = this->HowManyFracElement();
    int nel = fmeshvec[1]->NElements();
    for (int iel = 0; iel < nel; iel++) {
        TPZCompEl * cel = fmeshvec[1]->Element(iel);
        if (!cel) continue;
        TPZGeoEl * gel = cel->Reference();
        
        if (gel->Dimension() != 1) {
            continue;
        }
        
        int side = 1;
        TPZGeoElSide gelside(gel,side);
        TPZGeoElSide neigh = gelside.Neighbour();
        
        // Seeking for condition 1
        while (gelside!= neigh) {
            if (neigh.Element()->Dimension()==0 && neigh.Element()->MaterialId() == 8) {
                break;
            }
            neigh = neigh.Neighbour();
        }
        if (gelside == neigh) {
            continue;
        }
        
        TPZBlock<STATE> & block = fmeshvec[1]->Block();
        TPZGeoElSide gelsideleft(gel,0);
        TPZGeoElSide neighTip = gelsideleft.Neighbour();
        
        // Seeking for condition 2
        while (gelsideleft != neighTip) {
            if (neighTip.Element()->Dimension() == 1 && neighTip.Element()->MaterialId() == 6 ) {
                break;
            }
            neighTip = neighTip.Neighbour();
        }
        if (gelsideleft == neighTip) {
            DebugStop();
        }
        TPZGeoEl * gelleft = neighTip.Element();
        TPZCompEl * celleft = gelleft->Reference();
        
        
        // Chanching value
#ifdef DEBUG
        if (celleft->NConnects() != 1 && cel->NConnects() != 1) {
            DebugStop();
        }
#endif
        TPZConnect &connectleft =  celleft->Connect(0);
        TPZConnect &connect =  cel->Connect(0);
        
        int seqleft = connectleft.SequenceNumber();
        int seq = connect.SequenceNumber();
#ifdef DEBUG
        if (block.Size(seqleft) != 1 && block.Size(seq) != 1) {
            DebugStop();
        }
#endif
        int posleft = block.Position(seqleft);
        int pos = block.Position(seq);
        fmeshvec[1]->Solution()(pos,0) = fmeshvec[1]->Solution()(posleft,0);
        
        if (nfracel < 5) {
            TPZBuildMultiphysicsMesh::TransferFromMeshes(fmeshvec, fcmeshMixed);
            an->LoadSolution(fcmeshMixed->Solution());
            return;
        }
        
        // Seeking for condition 3
        TPZGeoElSide gelsidesecondleft(gelleft,0);
        TPZGeoElSide neighsecondeleft = gelsidesecondleft.Neighbour();
        while (gelsidesecondleft != neighsecondeleft) {
            if (neighsecondeleft.Element()->Dimension() == 1 && neighsecondeleft.Element()->MaterialId() == 6 ) {
                break;
            }
            neighsecondeleft = neighsecondeleft.Neighbour();
        }
        if (gelsidesecondleft == neighsecondeleft) {
            DebugStop();
        }
        TPZGeoEl * gelsecondleft = neighsecondeleft.Element();
        TPZCompEl * celsecondleft = gelsecondleft->Reference();
#ifdef DEBUG
        if (celsecondleft->NConnects() != 1) {
            DebugStop();
        }
#endif
        TPZConnect &connectsecondleft = celsecondleft->Connect(0);
        int seqsecondleft = connectsecondleft.SequenceNumber();
        
#ifdef DEBUG
        if (block.Size(seqsecondleft) != 1) {
            DebugStop();
        }
#endif
        int possecondleft = block.Position(seqsecondleft);
        
        fmeshvec[1]->Solution()(posleft,0) = fmeshvec[1]->Solution()(possecondleft,0);
        
        TPZBuildMultiphysicsMesh::TransferFromMeshes(fmeshvec, fcmeshMixed);
        an->LoadSolution(fcmeshMixed->Solution());
        return;
    }
}

int TPZDarcyAnalysis::HowManyFracElement()
{
    long nel = fgmesh->NElements();
    long nfracel = 0;
    for (long iel = 0; iel < nel; iel++) {
        TPZGeoEl *gel = fgmesh->Element(iel);
        if (gel->MaterialId() == 6) {
            nfracel++;
        }
    }
    return nfracel;
}

void TPZDarcyAnalysis::CreateMultiphysicsMesh(TPZFMatrix<REAL> Vl)
{
    fmeshvec[0] = CreateCMeshFluxHdiv();
    fmeshvec[1] = CreateCMeshPressureL2();
    
    //    // Initial Pressure
    //    TPZVec<STATE> solini(1,0.0);
    //    TPZCompMesh  * cmeshL2 = L2ProjectionP(fgmesh, fData->PorderPressure(), solini);
    //    TPZAnalysis anL2(cmeshL2);
    //    SolveSyst(anL2, cmeshL2);
    //    fmeshvec[1]->LoadSolution(anL2.Solution());
    
    fcmeshMixed = CreateCMeshMixed(Vl);
    
    // Transferindo para a multifisica
    TPZBuildMultiphysicsMesh::AddElements(fmeshvec, fcmeshMixed);
    TPZBuildMultiphysicsMesh::AddConnects(fmeshvec, fcmeshMixed);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(fmeshvec, fcmeshMixed);
    
    // Create Interfaces
    CreateInterfaces(fcmeshMixed);
    
    // Preparando os index dos pontos de integracao.
    long nel = fcmeshMixed->NElements();
    for (long iel = 0; iel < nel; iel++) {
        TPZCompEl *cel = fcmeshMixed->ElementVec()[iel];
        cel->PrepareIntPtIndices();
    }
    
    
    std::ofstream dumpfile("ComputationaMeshMultiphysics.txt");
    fcmeshMixed->Print(dumpfile);
}

void TPZDarcyAnalysis::AcceptSolution(TPZAnalysis *an)
{
    //Here mat 6 means fracture elements
    TPZMaterial *Material = fcmeshMixed->FindMaterial(6);
    TPZMatfrac1dhdiv * MaterialOfFract = dynamic_cast<TPZMatfrac1dhdiv *>(Material);
    MaterialOfFract->SetUpdateMem();
    an->AssembleResidual();
    MaterialOfFract->SetUpdateMem(false);
}

void TPZDarcyAnalysis::ComputeFirstSolForOneELement(TPZAnalysis * an)
{
    long nfracel = this->HowManyFracElement();
    
    if (nfracel != 1) {
        PZError << "This method sould only be called when the mesh has a single frac element " << std::endl;
        DebugStop();
    }
    
    int nel = fgmesh->NElements();
    TPZGeoEl *gel = NULL;
    for (int iel = 0; iel < nel; iel++) {
        gel = fgmesh->Element(iel);
        if (gel->MaterialId() == 6) {
            break;
        }
    }
    
    TPZBlock<STATE> &blockQ = fmeshvec[0]->Block(), &blockP = fmeshvec[1]->Block();
    
    
    // Setando os valores dos fluxos na fratura
    fgmesh->ResetReference();
    fmeshvec[0]->LoadReferences();
    
    const REAL pfrac = fData->SigmaConf();
    const REAL tstar = fData->FictitiousTime(fData->AccumVl(), pfrac); // VlForNextPropag is vl from last propag here
    const REAL vlnext = fData->VlFtau(pfrac, tstar+fData->TimeStep());
    const REAL totalLeakOff = 2. * fData->ElSize() * vlnext;
    const REAL totalLeakOffPrev = 2. * fData->ElSize() * fData->AccumVl();
    const REAL ql = (totalLeakOff - totalLeakOffPrev)/fData->TimeStep();
    REAL qout = fData->Q() - ql;
    
    TPZCompEl *celQ = gel->Reference();
    if (celQ->NConnects() != 3) {
        DebugStop(); // Mesh H1 1D p = 1
    }
    TPZConnect &c1Q = celQ->Connect(0), &c2Q = celQ->Connect(1);
    int seq1Q = c1Q.SequenceNumber(), seq2Q = c2Q.SequenceNumber();
    int pos1Q = blockQ.Position(seq1Q), pos2Q = blockQ.Position(seq2Q);
    fmeshvec[0]->Solution()(pos1Q,0) = fData->Q();
    fmeshvec[0]->Solution()(pos2Q,0) = qout;
    
    
    // Setando a pressao
    fgmesh->ResetReference();
    fmeshvec[1]->LoadReferences();
    const REAL dwdp = fData->GetDwDp();
    const REAL pini = fData->SigmaConf() + pow(12. * fData->Viscosity() * qout * fData->ElSize() / (dwdp*dwdp*dwdp),1./4.);
    
    TPZCompEl *celP = gel->Reference();
    if (celP->NConnects() != 1) {
        DebugStop(); // Mesh L2 1D p = 0
    }
    TPZConnect &cP = celP->Connect(0);
    int seqP = cP.SequenceNumber();
    int posP = blockP.Position(seqP);
    
    fmeshvec[1]->Solution()(posP,0) = pini;
    
    TPZBuildMultiphysicsMesh::TransferFromMeshes(fmeshvec, fcmeshMixed);
    an->LoadSolution(fcmeshMixed->Solution());
    
}

void TPZDarcyAnalysis::SetPressureOnNewElement(TPZAnalysis *an)
{
    fgmesh->ResetReference();
    fmeshvec[1]->LoadReferences();
    
    long nfracel = this->HowManyFracElement();
    int nel = fmeshvec[1]->NElements();
    for (int iel = 0; iel < nel; iel++) {
        TPZCompEl * cel = fmeshvec[1]->Element(iel);
        if (!cel) continue;
        TPZGeoEl * gel = cel->Reference();
        
        if (gel->Dimension() != 1) {
            continue;
        }
        
        int side = 1;
        TPZGeoElSide gelside(gel,side);
        TPZGeoElSide neigh = gelside.Neighbour();
        
        // Seeking for condition 1
        while (gelside!= neigh) {
            if (neigh.Element()->Dimension()==0 && neigh.Element()->MaterialId() == 8) {
                break;
            }
            neigh = neigh.Neighbour();
        }
        if (gelside == neigh) {
            continue;
        }
        
        TPZBlock<STATE> & block = fmeshvec[1]->Block();
        TPZGeoElSide gelsideleft(gel,0);
        TPZGeoElSide neighTip = gelsideleft.Neighbour();

#ifdef DEBUG
        if (cel->NConnects() != 1) {
            DebugStop();
        }
#endif
        TPZConnect &connect =  cel->Connect(0);
        int seq = connect.SequenceNumber();
#ifdef DEBUG
        if (block.Size(seq) != 1) {
            DebugStop();
        }
#endif
        int pos = block.Position(seq);
        fmeshvec[1]->Solution()(pos,0) = fData->SigmaConf();
        
        TPZBuildMultiphysicsMesh::TransferFromMeshes(fmeshvec, fcmeshMixed);
        an->LoadSolution(fcmeshMixed->Solution());
        return;
    }
}
