#include "Common.h"

#include "pzskylstrmatrix.h"
#include "TPZSSpStructMatrix.h"
#include "pzstepsolver.h"

#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif

#include "TPZMatElasticity2D.h"
#include "TPZMatLaplacian.h"
#include "pzbndcond.h"

#include "pzgengrid.h"
#include "TPZBuildSBFem.h"

#include "TPZVTKGeoMesh.h"

#include "TPZSBFemElementGroup.h"
#include "pzinterpolationspace.h"

#include "tpzarc3d.h"
#include "tpzgeoblend.h"

#include "pzcondensedcompel.h"

#ifdef _AUTODIFF
TElasticity2DAnalytic ElastExact;

TLaplaceExample1 LaplaceExact;

#endif

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.sbfem"));
//static LoggerPtr loggerBF(Logger::getLogger("pz.mesh.sbfemelementgroupBF"));
#endif

//TElasticity2DAnalytic::EDefState TElasticity2DAnalytic::fProblemType = TElasticity2DAnalytic::EStretchx;

int fExample;

void SolveSist(TPZAnalysis *an, TPZCompMesh *Cmesh)
{
    //    TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(Cmesh);
    TPZSkylineStructMatrix strmat(Cmesh);
    //    TPZSymetricSpStructMatrix strmat(Cmesh);
    strmat.SetNumThreads(0);
    an->SetStructuralMatrix(strmat);
    
    int64_t neq = Cmesh->NEquations();
    
    if(neq > 20000)
    {
        std::cout << "Entering Assemble Equations\n";
        std::cout.flush();
    }
#ifdef USING_BOOST
    boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
#endif
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an->SetSolver(step);
    
    an->Assemble();
    
    //    std::ofstream andrade("../Andrade.mtx");
    //    andrade.precision(16);
    //    an->Solver().Matrix()->Print("Andrade",andrade,EMatrixMarket);
    //    std::cout << "Leaving Assemble\n";
#ifdef USING_BOOST
    boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
#endif
    
    
    if(neq > 20000)
    {
        std::cout << "Entering Solve\n";
        std::cout.flush();
    }
    an->Solve();
    
#ifdef USING_BOOST
    boost::posix_time::ptime t3 = boost::posix_time::microsec_clock::local_time();
    std::cout << "Time for assembly " << t2-t1 << " Time for solving " << t3-t2 << std::endl;
#endif
    
    
}

//void Constant(const TPZVec<REAL> &x, TPZVec<STATE> &val)
//{
//    val.resize(1);
//    val[0] = 1.;
//}

void HarmonicNeumannLeft(const TPZVec<REAL> &x, TPZVec<STATE> &val)
{
    val[0] = -M_PI*exp(M_PI*x[0])*sin(M_PI*x[1]);
}

void HarmonicNeumannRight(const TPZVec<REAL> &x, TPZVec<STATE> &val)
{
    val[0] = M_PI*exp(M_PI*x[0])*sin(M_PI*x[1]);
}

void Harmonic_exact(const TPZVec<REAL> &xv, TPZVec<STATE> &val, TPZFMatrix<STATE> &deriv)
{
    val[0] = exp(M_PI*xv[0])*sin(M_PI*xv[1]);
    deriv(0,0) = M_PI*val[0];
    deriv(1,0) = M_PI*exp(M_PI*xv[0])*cos(M_PI*xv[1]);
    
}

void InsertMaterialObjects(TPZCompMesh *cmesh, bool scalarproblem, bool applyexact)
{
    
    // Getting mesh dimension
    int dim = 2;
    
    TPZMaterial *material;
    int nstate = 1;
    bool elasticity = false;
    if (!scalarproblem) {
        elasticity = true;
    }
    if (elasticity)
    {
        TPZMatElasticity2D *matloc1 = new TPZMatElasticity2D(Emat1);
        TPZMatElasticity2D *matloc2 = new TPZMatElasticity2D(Emat2);
        material = matloc1;
        nstate = 2;
        // Plane strain assumption
        //        REAL lamelambda = 1.0e9,lamemu = 0.5e3, fx= 0, fy = 0;
//        REAL lamelambda = 0.,lamemu = 0.5e3, fx= 0, fy = -0.5e3;
//        matloc1->SetParameters(lamelambda,lamemu, fx, fy);
//        matloc2->SetParameters(lamelambda,lamemu, fx, fy);
        TPZManVector<REAL,3> x(3,0.);
#ifdef _AUTODIFF
        // Setting up paremeters
        if (applyexact)
        {
            matloc1->SetPlaneStress();
            matloc2->SetPlaneStress();
//            matloc1->SetPlaneStrain();
//            matloc2->SetPlaneStrain();
            
            matloc1->SetElasticParameters(ElastExact.fE, ElastExact.fPoisson);
            matloc2->SetElasticParameters(ElastExact.fE, ElastExact.fPoisson);
        }
#endif
//        REAL Sigmaxx = 0.0, Sigmayx = 0.0, Sigmayy = 0.0, Sigmazz = 0.0;
//        matloc1->SetPreStress(Sigmaxx,Sigmayx,Sigmayy,Sigmazz);
//        matloc2->SetPreStress(Sigmaxx,Sigmayx,Sigmayy,Sigmazz);
        
#ifdef _AUTODIFF
        if(applyexact)
        {
            matloc1->SetForcingFunction(ElastExact.ForcingFunction());
            matloc2->SetForcingFunction(ElastExact.ForcingFunction());
        }
#endif
    }
    else
    {
        TPZMatLaplacian *matloc = new TPZMatLaplacian(Emat1);
        matloc->SetDimension(2);
        matloc->SetSymmetric();
        material = matloc;
        nstate = 1;
    }
    
    
    TPZDummyFunction<STATE> *dummy;
    TPZAutoPointer<TPZFunction<STATE> > autodummy;
    dummy = dummy = new TPZDummyFunction<STATE>(Elasticity_exact, 6);
    autodummy = dummy;
    
    TPZFMatrix<STATE> val1(nstate,nstate,0.), val2(nstate,1,0.);
//    val2(0,0) = 0.0;
//    val2(1,0) = 0.0;
    TPZMaterial * BCond1;
    if(!elasticity)
    {
        BCond1 = material->CreateBC(material,Ebc1,0, val1, val2);
    }
    else
    {
//        BCond1 = material->CreateBC(material,Ebc1,1, val1, val2);
        BCond1 = material->CreateBC(material,Ebc1,0, val1, val2);
#ifdef _AUTODIFF
        if (applyexact) {
//            BCond1->SetForcingFunction(ElastExact.TensorFunction());
            BCond1->SetForcingFunction(autodummy);
        }
#endif
    }
    
//    val2(0,0) = -1.0*1000.0;
//    val2(1,0) = 0.0;
    TPZMaterial * BCond2 = 0;
    if(elasticity == 0)
    {
        BCond2 = material->CreateBC(material,Ebc2,0, val1, val2);
    }
    else
    {
        // mixed condition on the right side to desingularize the problem
//        val1(0,0) = 1.;
//        val1(1,1) = 1.;
//        BCond2 = material->CreateBC(material,Ebc2,1, val1, val2);
        BCond2 = material->CreateBC(material,Ebc2,0, val1, val2);
#ifdef _AUTODIFF
        if (applyexact) {
//            BCond2->SetForcingFunction(ElastExact.TensorFunction());
            BCond2->SetForcingFunction(autodummy);
        }
#endif
        val1.Zero();
    }
//    val2(0,0) = 0.0;
//    val2(1,0) = 0.0;
    TPZMaterial * BCond3;
    if(!elasticity)
    {
        BCond3 = material->CreateBC(material,Ebc3,0, val1, val2);
    }
    else
    {
//        BCond3 = material->CreateBC(material,Ebc3,1, val1, val2);
        BCond3 = material->CreateBC(material,Ebc3,0, val1, val2);
#ifdef _AUTODIFF
        if (applyexact) {
//            BCond3->SetForcingFunction(ElastExact.TensorFunction());
            BCond3->SetForcingFunction(autodummy);
        }
#endif
    }
    
//    val2(0,0) = -1.0*1000.0;
//    if(elasticity) val2(0,0) *=-1.;
//    TPZMaterial * BCond4 = material->CreateBC(material,Ebc4,1, val1, val2);
    TPZMaterial * BCond4 = material->CreateBC(material,Ebc4,0, val1, val2);
#ifdef _AUTODIFF
    if (applyexact) {
//        BCond4->SetForcingFunction(ElastExact.TensorFunction());
        BCond4->SetForcingFunction(autodummy);
    }
    
//    if (elasticity && applyexact) {
//        val1.Zero();
//        val1(0,0) = 0.01;
//        val1(1,1) = 0.01;
//        TPZMaterial * BCond5 = material->CreateBC(material,EBCPoint1, 0, val1, val2);
////        BCond5->SetForcingFunction(ElastExact.Exact());
//        val1(0,0) = 0.;
//        TPZMaterial * BCond6 = material->CreateBC(material,EBCPoint2, 0, val1, val2);
////        BCond6->SetForcingFunction(ElastExact.Exact());
//        BCond5->SetForcingFunction(autodummy);
//        BCond5->SetForcingFunction(autodummy);
//        cmesh->InsertMaterialObject(BCond5);
//        cmesh->InsertMaterialObject(BCond6);
//    }
#endif
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZMaterial * BSkeleton = material->CreateBC(material,ESkeleton,1, val1, val2);
    
    
    cmesh->InsertMaterialObject(material);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
    cmesh->InsertMaterialObject(BSkeleton);
    
}



void InsertMaterialObjects2(TPZCompMesh *cmesh, bool scalarproblem, bool applyexact)
{
    
    // Getting mesh dimension
    int dim = 2;
    
    TPZMaterial *material;
    int nstate = 1;
    bool elasticity = false;
    if (!scalarproblem) {
        elasticity = true;
    }
    
    int polynomialorder = 6;
    TPZDummyFunction<STATE> *dummy;
    TPZDummyFunction<STATE> *dummyforce;
    TPZAutoPointer<TPZFunction<STATE> > autodummy;
    TPZAutoPointer<TPZFunction<STATE> > autodummyforce;
    if (!elasticity){
        switch (fExample) {
            case 1:
                dummy = new TPZDummyFunction<STATE>(ExactSol_Laplacian_Ex1, polynomialorder);
                break;
                
            case 2:
                dummy = new TPZDummyFunction<STATE>(ExactSol_Laplacian_Ex2, polynomialorder);
                break;
                
            case 3:
                dummy = new TPZDummyFunction<STATE>(ExactSol_Laplacian_Ex3, polynomialorder);
                break;
                
            case 4:
                dummy = new TPZDummyFunction<STATE>(ExactSol_Laplacian_Ex4, polynomialorder);
                break;
                
            case 5:
                dummy = new TPZDummyFunction<STATE>(Laplace_exact, polynomialorder);
                break;
                
            default:
                DebugStop();
                break;
        }
        autodummy = dummy;
      
        switch (fExample) {
            case 1:
                dummyforce = new TPZDummyFunction<STATE>(Constant_Laplacian_Ex1, polynomialorder);
                break;
                
            case 2:
                dummyforce = new TPZDummyFunction<STATE>(Constant_Laplacian_Ex2, polynomialorder);
                break;
                
            case 3:
                dummyforce = new TPZDummyFunction<STATE>(Constant_Laplacian_Ex3, polynomialorder);
                break;
                
            case 4:
                dummyforce = new TPZDummyFunction<STATE>(Constant_Laplacian_Ex4, polynomialorder);
                break;
                
            case 5:
                dummyforce = new TPZDummyFunction<STATE>(Laplace_exact, polynomialorder);
                break;
                
            default:
                DebugStop();
                break;
        }
        autodummyforce = dummyforce;
        if (fExample == 5) {
            autodummyforce = LaplaceExact.ForcingFunction();
        }
    }
    
    
    if (elasticity)
    {
        TPZMatElasticity2D *matloc1 = new TPZMatElasticity2D(Emat1);
        TPZMatElasticity2D *matloc2 = new TPZMatElasticity2D(Emat2);
        material = matloc1;
        nstate = 2;
        // Plane strain assumption
//                REAL lamelambda = 1.0e9,lamemu = 0.5e3, fx= 0, fy = 0;
        REAL lamelambda = 0.,lamemu = 0.5e3, fx= 0, fy = 0;
        matloc1->SetParameters(lamelambda,lamemu, fx, fy);
        matloc2->SetParameters(lamelambda,lamemu, fx, fy);
        TPZManVector<REAL,3> x(3,0.);
#ifdef _AUTODIFF
        // Setting up paremeters
        if (applyexact)
        {
            matloc1->SetPlaneStress();
            matloc1->SetElasticParameters(ElastExact.fE, ElastExact.fPoisson);
            matloc2->SetPlaneStress();
            matloc2->SetElasticParameters(ElastExact.fE, ElastExact.fPoisson);
        }
#endif
        REAL Sigmaxx = 0.0, Sigmayx = 0.0, Sigmayy = 0.0, Sigmazz = 0.0;
        matloc1->SetPreStress(Sigmaxx,Sigmayx,Sigmayy,Sigmazz);
        matloc2->SetPreStress(Sigmaxx,Sigmayx,Sigmayy,Sigmazz);
        
#ifdef _AUTODIFF
        if(applyexact)
        {
            matloc1->SetForcingFunction(ElastExact.ForcingFunction());
            matloc2->SetForcingFunction(ElastExact.ForcingFunction());
        }
#endif
    }
    else
    {
        TPZMatLaplacian *matloc = new TPZMatLaplacian(Emat1);
        matloc->SetDimension(2);
        matloc->SetSymmetric();
        material = matloc;
        nstate = 1;
        
        matloc->SetForcingFunction(autodummyforce);
    }
    
    
    
    TPZFMatrix<STATE> val1(nstate,nstate,0.), val2(nstate,1,0.);
    TPZMaterial * BCond1;
    if(!elasticity)
    {
        BCond1 = material->CreateBC(material,Ebc1,0, val1, val2);
        BCond1->SetForcingFunction(autodummy);
    }
    else
    {
        BCond1 = material->CreateBC(material,Ebc1,1, val1, val2);
#ifdef _AUTODIFF
        if (applyexact) {
            BCond1->SetForcingFunction(ElastExact.TensorFunction());
        }
#endif
    }
    
    val1.Zero();
    val2.Zero();
    val2(0,0) = -1.0*1000.0;
    TPZMaterial * BCond2 = 0;
    if(elasticity == 0)
    {
        BCond2 = material->CreateBC(material,Ebc2,0, val1, val2);
        BCond2->SetForcingFunction(autodummy);
    }
    else
    {
        // mixed condition on the right side to desingularize the problem
//        val1(0,0) = 1.;
//        val1(1,1) = 1.;
        BCond2 = material->CreateBC(material,Ebc2,1, val1, val2);
#ifdef _AUTODIFF
        if (applyexact) {
            BCond2->SetForcingFunction(ElastExact.TensorFunction());
        }
#endif
        val1.Zero();
    }
    
    
    val1.Zero();
    val2.Zero();
    TPZMaterial * BCond3;
    if(!elasticity)
    {
        BCond3 = material->CreateBC(material,Ebc3,0, val1, val2);
        BCond3->SetForcingFunction(autodummy);
    }
    else
    {
        BCond3 = material->CreateBC(material,Ebc3,1, val1, val2);
#ifdef _AUTODIFF
        if (applyexact) {
            BCond3->SetForcingFunction(ElastExact.TensorFunction());
        }
#endif
    }
    
    val1.Zero();
    val2.Zero();
    val2(0,0) = -1.0*1000.0;
    TPZMaterial * BCond4;
    if(!elasticity)
    {
        BCond4 = material->CreateBC(material,Ebc4,0, val1, val2);
        BCond4->SetForcingFunction(autodummy);
    }
    else
    {
        BCond4 = material->CreateBC(material,Ebc4,1, val1, val2);
#ifdef _AUTODIFF
        if (applyexact) {
            BCond4->SetForcingFunction(ElastExact.TensorFunction());
        }
#endif
    }
    
    if (elasticity && applyexact) {
        val1.Zero();
        val1(0,0) = 0.01;
        val1(1,1) = 0.01;
        TPZMaterial * BCond5 = material->CreateBC(material,EBCPoint1, 0, val1, val2);
        BCond5->SetForcingFunction(ElastExact.TensorFunction());
        val1(0,0) = 0.;
        TPZMaterial * BCond6 = material->CreateBC(material,EBCPoint2, 0, val1, val2);
        BCond6->SetForcingFunction(ElastExact.TensorFunction());
        cmesh->InsertMaterialObject(BCond5);
        cmesh->InsertMaterialObject(BCond6);
    }
    
    val1.Zero(); val2.Zero();
    TPZMaterial * BSkeleton = material->CreateBC(material,ESkeleton,1, val1, val2);
    
    cmesh->InsertMaterialObject(material);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
    cmesh->InsertMaterialObject(BSkeleton);
    
}



TPZCompMesh *SetupSquareMesh(int nelx, int nrefskeleton, int porder, bool scalarproblem, bool useexact)
{
    bool elasticityproblem = !scalarproblem;
    TPZManVector<REAL,4> x0(3,-1.),x1(3,1.);
//    x0[0] = -1;
//    x0[1] = -1;
    x0[0] = 0;
    x0[1] = 0;
    x1[0] = 1;
    x1[1] = 1;
    x0[2] = 0.;
    x1[2] = 0.;
    TPZManVector<int,4> nx(2,nelx);
    TPZGenGrid gengrid(nx,x0,x1);
    gengrid.SetElementType(EQuadrilateral);
    TPZAutoPointer<TPZGeoMesh> gmesh = new TPZGeoMesh;
    
    //        OneQuad(gmesh);
    gengrid.Read(gmesh,EGroup);
    gengrid.SetBC(gmesh, 4, Ebc1);
    gengrid.SetBC(gmesh, 5, Ebc2);
    gengrid.SetBC(gmesh, 6, Ebc3);
    gengrid.SetBC(gmesh, 7, Ebc4);
    {
        TPZManVector<int64_t,2> nodeindex(1);
        int64_t index;
        nodeindex[0] = 0;
        gmesh->CreateGeoElement(EPoint, nodeindex, EBCPoint1, index);
        nodeindex[0] = nelx;
        gmesh->CreateGeoElement(EPoint, nodeindex, EBCPoint2, index);
        gmesh->BuildConnectivity();
    }
    
    std::map<int,int> matmap;
    matmap[EGroup] = 1;
    TPZBuildSBFem build(gmesh,ESkeleton,matmap);
    
    build.StandardConfiguration();
    build.DivideSkeleton(nrefskeleton);
    //        AddSkeletonElements(gmesh);
    /// generate the SBFem elementgroups
    
    /// put sbfem pyramids into the element groups
    TPZCompMesh *SBFem = new TPZCompMesh(gmesh);
    SBFem->SetDefaultOrder(porder);
    
    // problemtype - 1 laplace equation
    int problemtype  = 0;
    if (elasticityproblem) {
        problemtype = 0;
    }
    else
    {
        problemtype = 1;
    }
    if(elasticityproblem) {
        InsertMaterialObjects(SBFem,!elasticityproblem, useexact);
    } else{
        InsertMaterialObjects2(SBFem,!elasticityproblem, useexact);
    }
//    if(problemtype == 1)
//    {
//        TPZMaterial *BCond2 = SBFem->FindMaterial(Ebc2);
//        TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(HarmonicNeumannRight,0);
//        TPZAutoPointer<TPZFunction<STATE> > autodummy = dummy;
//        BCond2->SetForcingFunction(autodummy);
//    }
//    if(problemtype == 1)
//    {
//        TPZMaterial *BCond4 = SBFem->FindMaterial(Ebc4);
//        TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(HarmonicNeumannLeft,0);
//        TPZAutoPointer<TPZFunction<STATE> > autodummy = dummy;
//        BCond4->SetForcingFunction(autodummy);
//    }
    
    
    build.BuildComputationMesh(*SBFem);
    
    if(1)
    {
        std::ofstream outc("CMesh.txt");
        SBFem->Print(outc);
        std::ofstream outg("GMesh.txt");
        gmesh->Print(outg);
        std::ofstream out("Geometry.vtk");
        TPZVTKGeoMesh vtk;
        vtk.PrintGMeshVTK(gmesh, out,true);
    }
    return SBFem;
}

TPZCompMesh *SetupSquareH1Mesh(int nelx, int porder, bool scalarproblem, bool useexact)
{
    bool elasticityproblem = !scalarproblem;
    TPZManVector<REAL,4> x0(3,-1.),x1(3,1.);
    x0[0] = -1;
    x0[1] = -1;
    x1[0] = 1;
    x1[1] = 1;
    x0[2] = 0.;
    x1[2] = 0.;
    TPZManVector<int,4> nx(2,nelx);
    TPZGenGrid gengrid(nx,x0,x1);
    gengrid.SetElementType(EQuadrilateral);
    TPZAutoPointer<TPZGeoMesh> gmesh = new TPZGeoMesh;
    
    //        OneQuad(gmesh);
    gengrid.Read(gmesh,Emat1);
    {
        TPZManVector<int64_t,2> nodeindex(1);
		int64_t index;
        nodeindex[0] = 0;
        gmesh->CreateGeoElement(EPoint, nodeindex, EBCPoint1, index);
        nodeindex[0] = nelx;
        gmesh->CreateGeoElement(EPoint, nodeindex, EBCPoint2, index);
    }
    gmesh->BuildConnectivity();
    
    gengrid.SetBC(gmesh, 4, Ebc1);
    gengrid.SetBC(gmesh, 5, Ebc2);
    gengrid.SetBC(gmesh, 6, Ebc3);
    gengrid.SetBC(gmesh, 7, Ebc4);
    gmesh->BuildConnectivity();
    
    /// put sbfem pyramids into the element groups
    TPZCompMesh *SBFem = new TPZCompMesh(gmesh);
    SBFem->SetDefaultOrder(porder);
    
    // problemtype - 1 laplace equation
    int problemtype  = 0;
    if (elasticityproblem) {
        problemtype = 0;
    }
    else
    {
        problemtype = 1;
    }
    
    if(elasticityproblem) {
        InsertMaterialObjects(SBFem,!elasticityproblem, useexact);
    } else{
        InsertMaterialObjects2(SBFem,!elasticityproblem, useexact);
    }
    
    SBFem->SetAllCreateFunctionsContinuous();
    SBFem->AutoBuild();
    
    if(1)
    {
        std::ofstream outc("CMesh.txt");
        SBFem->Print(outc);
        std::ofstream outg("GMesh.txt");
        gmesh->Print(outg);
        std::ofstream out("Geometry.vtk");
        TPZVTKGeoMesh vtk;
        vtk.PrintGMeshVTK(gmesh, out,true);
    }
    return SBFem;
}


TPZCompMesh *SetupOneArc(int numrefskeleton, int porder, REAL angle)
{
    TPZAutoPointer<TPZGeoMesh> gmesh = new TPZGeoMesh;
    gmesh->SetDimension(2);
    gmesh->NodeVec().Resize(4);
    TPZManVector<REAL,3> co(3,0.);
    gmesh->NodeVec()[0].Initialize(co, gmesh);
    co[0] = 1.;
    gmesh->NodeVec()[1].Initialize(co, gmesh);
    co.Fill(0.);
    co[0] = cos(angle);
    co[1] = sin(angle);
    gmesh->NodeVec()[2].Initialize(co, gmesh);
    co.Fill(0.);
    co[0] = cos(angle/2.);
    co[1] = sin(angle/2.);
    gmesh->NodeVec()[3].Initialize(co, gmesh);
    co.Fill(0.);
    TPZManVector<int64_t,4> nodeindex(1,0);
    
    nodeindex[0] = 1;
    int64_t elementid = 1;
    gmesh->CreateGeoElement(EPoint, nodeindex, Ebc1, elementid);
    
    nodeindex.Resize(3);
    // Definition of Arc coordenates
    // Create Geometrical Arc #1
    nodeindex[0] = 1;
    nodeindex[1] = 2;
    nodeindex[2] = 3;
    elementid = 1;
    TPZGeoEl *arc = new TPZGeoElRefPattern < pzgeom::TPZArc3D > (nodeindex, Ebc2, gmesh,elementid);
    
    nodeindex.Resize(4);
    nodeindex[0] = 1;
    nodeindex[1] = 2;
    nodeindex[2] = 0;
    nodeindex[3] = 0;
    elementid = 2;
    TPZGeoEl *gblend = new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> > (nodeindex, EGroup, gmesh,elementid);
    
    gmesh->BuildConnectivity();
    
    //    gmesh->Print(std::cout);
    TPZManVector<REAL,3> xi(1),x(3);
    for (REAL s=-1.; s<=1.; s+= 1./10.) {
        xi[0] = s;
        arc->X(xi, x);
        std::cout << "xi " << xi << " x " << x << std::endl;
    }
    std::map<int,int> matmap;
    matmap[EGroup] = Emat1;
    TPZBuildSBFem build(gmesh,ESkeleton,matmap);
    
    TPZManVector<int64_t,5> elids(1,gblend->Index());
    build.AddPartition(elids, 0);
    
    build.DivideSkeleton(numrefskeleton);
    //        AddSkeletonElements(gmesh);
    /// generate the SBFem elementgroups
    
    /// put sbfem pyramids into the element groups
    TPZCompMesh *SBFem = new TPZCompMesh(gmesh);
    SBFem->SetDefaultOrder(porder);
    
    // problemtype - 1 laplace equation
    int problemtype  = 1;
    bool applyexact = false;
    InsertMaterialObjects(SBFem,problemtype,applyexact);
    
    
    build.BuildComputationMesh(*SBFem);
    
    {
        std::ofstream outg("GMesh.txt");
        gmesh->Print(outg);
        std::ofstream outc("CMesh.txt");
        SBFem->Print(outc);
        std::ofstream out("Geometry.vtk");
        TPZVTKGeoMesh vtk;
        vtk.PrintGMeshVTK(gmesh, out,true);
    }
    return SBFem;
    
}

void ElGroupEquations(TPZSBFemElementGroup *elgr, TPZVec<int64_t> &equations)
{
    equations.Resize(0, 0);
    TPZCompMesh *cmesh = elgr->Mesh();
    int nc = elgr->NConnects();
    for (int ic=0; ic<nc; ic++) {
        TPZConnect &c = elgr->Connect(ic);
        int blsize = c.NDof();
        int64_t eqsize = equations.size();
        equations.Resize(eqsize+blsize, 0);
        int64_t seqnum = c.SequenceNumber();
        for (int idf = 0; idf<blsize; idf++) {
            equations[eqsize+idf] = cmesh->Block().Position(seqnum)+idf;
        }
    }
}
/// Verify if the values of the shapefunctions corresponds to the value of ComputeSolution for all SBFemVolumeElements
void VerifyShapeFunctionIntegrity(TPZSBFemVolume *celv)
{
    TPZGeoEl *gel = celv->Reference();
    int dim = gel->Dimension();
    int nstate = celv->Connect(0).NState();
    TPZCompMesh *cmesh = celv->Mesh();
    int volside = gel->NSides()-1;
    TPZSBFemElementGroup *elgr = dynamic_cast<TPZSBFemElementGroup *>(cmesh->Element(celv->ElementGroupIndex()));
    TPZManVector<int64_t> globeq;
    ElGroupEquations(elgr, globeq);
    TPZIntPoints *intpoints = gel->CreateSideIntegrationRule(volside, 3);
    cmesh->Solution().Zero();
    for (int ip=0; ip < intpoints->NPoints(); ip++) {
        TPZManVector<REAL,3> xi(gel->Dimension(),0.);
        TPZFNMatrix<32,REAL> phi,dphidxi;
        REAL weight;
        intpoints->Point(ip, xi, weight);
        celv->Shape(xi, phi, dphidxi);
        int64_t neq = globeq.size();
        for (int64_t eq=0; eq<neq; eq++) {
            int64_t globindex = globeq[eq];
            cmesh->Solution().Zero();
            cmesh->Solution()(globindex,0) = 1.;
            cmesh->LoadSolution(cmesh->Solution());
            TPZSolVec sol;
            TPZGradSolVec dsol;
            TPZFNMatrix<9,REAL> axes(dim,3);
            celv->ComputeSolution(xi, sol, dsol, axes);
            REAL diffphi = 0., diffdphi = 0.;
            for (int istate = 0; istate < nstate; istate++) {
                diffphi += (sol[0][istate]-phi(eq*nstate+istate))*(sol[0][istate]-phi(eq*nstate+istate));
                for (int d=0; d<dim; d++) {
                    STATE diff = (dsol[0](d,istate)-dphidxi(d+istate*nstate,eq));
                    diffdphi += diff*diff;
                }
            }
            diffphi = sqrt(diffphi);
            diffdphi = sqrt(diffdphi);
            if (diffphi > 1.e-8 || diffdphi > 1.e-8) {
                std::cout << "Wrong shape function diffphi = " << diffphi << " diffdphi " << diffdphi << "\n";
            }
        }
    }
    delete intpoints;
}

/// Verify if the values of the shapefunctions corresponds to the value of ComputeSolution for all SBFemVolumeElements
void VerifyShapeFunctionIntegrity(TPZCompMesh *cmesh)
{
    int64_t nel = cmesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZSBFemElementGroup *elgr = dynamic_cast<TPZSBFemElementGroup *>(cel);
        if (elgr) {
            TPZStack<TPZCompEl *,5> elstack = elgr->GetElGroup();
            int nvol = elstack.size();
            for (int iv=0; iv<nvol; iv++) {
                TPZCompEl *vcel = elstack[iv];
                TPZSBFemVolume *elvol = dynamic_cast<TPZSBFemVolume *>(vcel);
                VerifyShapeFunctionIntegrity(elvol);
            }
        }
    }
}

/// Build a square mesh with boundary conditions
TPZCompMesh *SetupCrackedOneElement(int nrefskeleton, int porder, bool applyexact, bool elastic)
{
    TPZAutoPointer<TPZGeoMesh> gmesh = new TPZGeoMesh;
    gmesh->SetDimension(2);
    REAL coor[][3] = {
        {0,0},
        {-1,0},
        {-1,-1},
        {1,-1},
        {1,1},
        {-1,1},
        {-1,0}
    };
    gmesh->NodeVec().Resize(7);
    for (int i=0; i<7; i++) {
        TPZManVector<REAL,3> co(3,0);
        co[0] = coor[i][0];
        co[1] = coor[i][1];
        gmesh->NodeVec()[i].Initialize(co, gmesh);
    }
    {
        
        TPZManVector<int64_t,2> nodeindices(2);
        nodeindices[0] = 1;
        nodeindices[1] = 2;
        int64_t index;
        gmesh->CreateGeoElement(EOned, nodeindices, Emat1, index);
        gmesh->CreateGeoElement(EOned, nodeindices, Ebc1, index);
        for (int i=1; i<4; i++) {
            nodeindices[0] = i+1;
            nodeindices[1] = i+2;
            gmesh->CreateGeoElement(EOned, nodeindices, Emat2, index);
            gmesh->CreateGeoElement(EOned, nodeindices, Ebc2, index);
        }
        nodeindices[0] = 5;
        nodeindices[1] = 6;
        gmesh->CreateGeoElement(EOned, nodeindices, Emat3, index);
        gmesh->CreateGeoElement(EOned, nodeindices, Ebc3, index);
    }
    gmesh->BuildConnectivity();
    std::map<int,int> matidtranslation;
    matidtranslation[Emat1] = Emat1;
    matidtranslation[Emat2] = Emat2;
    matidtranslation[Emat3] = Emat3;
    TPZBuildSBFem build(gmesh, ESkeleton, matidtranslation);
    TPZManVector<int64_t,10> scalingcenters(1);
    scalingcenters[0] = 0;
    int64_t nel = gmesh->NElements();
    TPZManVector<int64_t,10> elementgroup(nel,-1);
    for (int i=0; i<nel; i+=2) {
        elementgroup[i] = 0;
    }
    build.SetPartitions(elementgroup, scalingcenters);
    std::set<int> matids;
    matids.insert(Ebc1);
    matids.insert(Ebc2);
    matids.insert(Ebc3);
    build.DivideSkeleton(nrefskeleton,matids);
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(porder);
    InsertMaterialObjects(cmesh, !elastic, true);
    TPZMaterial *mat = cmesh->FindMaterial(Emat1);
    TPZMaterial *mat2 = mat->NewMaterial();
    mat2->SetId(Emat2);
    cmesh->InsertMaterialObject(mat2);
    TPZMaterial *mat3 = mat->NewMaterial();
    mat3->SetId(Emat3);
    cmesh->InsertMaterialObject(mat3);
    build.BuildComputationalMeshFromSkeleton(*cmesh);
    {
        int64_t nel = gmesh->NElements();
        for (int64_t el=0; el<nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if (gel->Dimension() != gmesh->Dimension()-1) {
                continue;
            }
            if (gel->MaterialId() == Emat1 || gel->MaterialId() == Emat2 || gel->MaterialId() == Emat3) {
                gel->SetMaterialId(ESkeleton);
            }
        }
    }
    return cmesh;
}


void Constant_Laplacian_Ex1(const TPZVec<REAL> &x, TPZVec<STATE> &val)
{
    val[0] = 0.0;
}

void ExactSol_Laplacian_Ex1(const TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix<STATE> &deriv)
{
    sol.resize(1);
    deriv.Resize(2,1);
    
    sol[0] = 1./2*(-x[0]*x[0]+x[1]*x[1]);
    deriv(0,0) = -x[0];
    deriv(1,0) = x[1];
}


void Constant_Laplacian_Ex2(const TPZVec<REAL> &x, TPZVec<STATE> &val)
{
    val[0] = 0.0;
}

void ExactSol_Laplacian_Ex2(const TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix<STATE> &deriv)
{
    sol.resize(1);
    deriv.Resize(2,1);
    
    sol[0] = 1./2 * ( sin(M_PI*x[0])*sinh(M_PI*x[1]) );
    deriv(0,0) = 1./2 * M_PI * cos(M_PI*x[0]) * sinh(M_PI*x[1]);
    deriv(1,0) = 1./2 * M_PI * cosh(M_PI*x[1]) * sin(M_PI*x[0]);
}

void Constant_Laplacian_Ex3(const TPZVec<REAL> &x, TPZVec<STATE> &val)
{
    val[0] = 1.0;
}

void ExactSol_Laplacian_Ex3(const TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix<STATE> &deriv)
{
    sol.resize(1);
    deriv.Resize(2,1);
    
    sol[0] = 1./2*sin(M_PI*x[0])*sinh(M_PI*x[1]) - x[0]*x[0]/2;
    deriv(0,0) = 1./2*M_PI*cos(M_PI*x[0])*sinh(M_PI*x[1]) - x[0];
    deriv(1,0) = 1./2*M_PI*cosh(M_PI*x[1])*sin(M_PI*x[0]);
}

void Constant_Laplacian_Ex4(const TPZVec<REAL> &x, TPZVec<STATE> &val)
{
    val[0] = 1.0;
}

void ExactSol_Laplacian_Ex4(const TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix<STATE> &deriv)
{
    
    sol.resize(1);
    deriv.Resize(2,1);
    
    sol[0] = 1./2*(1.-x[0]*x[0]);
    deriv(0,0) = -x[0];
    deriv(1,0) = 0;
}



void BodyLoads_ElasticityOoietal1(const TPZVec<REAL> &cord, TPZVec<REAL> &val)
{
    val.resize(2);
    
    REAL x = cord[0];
    REAL y = cord[1];
    
    val[0] = (7*pow(y,2)*pow(-1 + 1.*y,3) + pow(x,5)*(2. - 18.*y + 36.*pow(y,2) - 20.*pow(y,3)) + x*y*(-10. + 108.*y - 249.*pow(y,2) + 214.*pow(y,3) - 63.*pow(y,4)) +
    pow(x,4)*(-6. + 79.*y - 220.5*pow(y,2) + 210.*pow(y,3) - 62.5*pow(y,4)) +
    pow(x,3)*(6. - 114.*y + 448.*pow(y,2) - 630.*pow(y,3) + 360.*pow(y,4) - 70.*pow(y,5)) +
    pow(x,2)*(-2. + 63.*y - 364.5*pow(y,2) + 668.*pow(y,3) - 490.5*pow(y,4) + 126.*pow(y,5)));
    
    val[1] = (2.*pow(y,2)*pow(-1. + 1.*y,3) + pow(x,5)*(7. - 63.*y + 126.*pow(y,2) - 70.*pow(y,3)) +
    pow(x,4)*(-21. + 214.*y - 490.5*pow(y,2) + 360.*pow(y,3) - 62.5*pow(y,4)) + x*y*(-10. + 63.*y - 114.*pow(y,2) + 79.*pow(y,3) - 18.*pow(y,4)) +
    pow(x,3)*(21. - 249.*y + 668.*pow(y,2) - 630.*pow(y,3) + 210.*pow(y,4) - 20.*pow(y,5)) +
    pow(x,2)*(-7. + 108.*y - 364.5*pow(y,2) + 448.*pow(y,3) - 220.5*pow(y,4) + 36.*pow(y,5)));
}

void ExactSol_ElasticityOoietal1(const TPZVec<REAL> &cord, TPZVec<REAL> &sol, TPZFMatrix<STATE> &deriv)
{
    
    sol.resize(2);
    deriv.Resize(2,2);
    
    REAL x = cord[0];
    REAL y = cord[1];
    
    sol[0] = x*x*y*y*pow(1-x,3)*pow(1-y,3);
    sol[1] = x*x*y*y*pow(1-x,3)*pow(1-y,3);
    
    deriv(0,0) = 2*pow(1-x,3)*x*pow(1-y,3)*y*y - 3*(1-x)*(1-x)*x*x*pow(1-y,3)*y*y;
    deriv(1,0) = 2*pow(1-x,3)*x*x*pow(1-y,3)*y - 3*pow(1-x,3)*x*x*(1-y)*(1-y)*y*y;
    deriv(0,1) = 2*pow(1-x,3)*x*pow(1-y,3)*y*y - 3*(1-x)*(1-x)*pow(1-y,3)*y*y;
    deriv(1,1) = 2*pow(1-x,3)*x*x*pow(1-y,3)*y - 3*pow(1-x,3)*x*x*(1-y)*(1-y)*y*y;
}

void forcefunction(const TPZVec<REAL> &co, TPZVec<STATE> &result)
{
    result.resize(1);
    
    double xv = co[0];
    double yv = co[1];
    
    double Pi = M_PI;
    result[0] = 8.0*Pi*Pi*cos(2.0*Pi*yv)*sin(2.0*Pi*xv);
}

void solexact(const TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix<STATE> &deriv)
{
    sol.resize(1);
    deriv.Resize(2,1);
    
    double xv = x[0];
    double yv = x[1];
    
    double Pi = M_PI;
    
    double v_x = cos(2*Pi*yv)*sin(2*Pi*xv);
    
    sol[0] = v_x;
//    result[1] = -v_x;
    
    deriv(0,0)=  2*Pi*cos(2*Pi*yv)*cos(2*Pi*xv);
    deriv(1,0)= -2*Pi*sin(2*Pi*xv)*sin(2*Pi*yv);
//    deriv(0,1)= -2*Pi*cos(2*Pi*yv)*cos(2*Pi*xv);
//    deriv(1,1)=  2*Pi*sin(2*Pi*xv)*sin(2*Pi*yv);
}

void WriteCSVFile(std::ofstream results, TPZVec<REAL> &errors)
{
    for(int i=0;i<3;i++) errors[i] *= 1e6;
    std::stringstream varname;
//    results << "Errmat[[" << nelxcount << ", " << POrder << ", " << 1 << ";;" << 3 << "]] = (1/1000000)*";\
//    results << "{ " << errors << " };" << std::endl;
    
}

void SetExample(int example){
    fExample = example;
}

void SetInterpolation(TPZCompMesh * cmesh, REAL h){
    
    int64_t nel = cmesh->NElements();

    for (int64_t iel=0; iel<nel; iel++) {
        if(!cmesh->Element(iel)){
            continue;
        }
        
        TPZCompEl * cel = cmesh->Element(iel);
        
        TPZSBFemElementGroup *celelem = dynamic_cast<TPZSBFemElementGroup *>(cel);
        if(!celelem){
            continue;
        }
        TPZStack<TPZCompEl *, 5> elgroup = celelem->GetElGroup();
        
        for (int i=0; i<elgroup.size(); i++) {
            TPZSBFemVolume *celvol = dynamic_cast<TPZSBFemVolume *>(elgroup[i]);
            if (!celvol) {
                continue;
            }
            TPZCompEl *celskeleton = cmesh->Element(celvol->SkeletonIndex());
            TPZGeoEl * gel = celskeleton->Reference();
            
//            TPZVec<REAL> qsi(3,0.5);
//            TPZFMatrix<REAL> jac;
//            TPZFMatrix<REAL> axes;
//            REAL detjac;
//            TPZFMatrix<REAL> jacinv;
//            gel->Jacobian(qsi, jac, axes, detjac, jacinv);
//            TPZFMatrix<REAL> coord;
//            gel->NodesCoordinates(coord);
//            
//            celvol->Coeficients().Print(std::cout);
//
//            if (coord(1,0) == coord(1,1)) {
//                TPZFNMatrix<200,std::complex<double>> coef(3,1,0);
//                TPZManVector<REAL> x(2,0);
//                TPZManVector<REAL> sol0(1,0);
//                TPZManVector<REAL> sol1(1,0);
//                TPZFMatrix<STATE> deriv;
//                
//                x[0] = coord(0,0);
//                x[1] = coord(1,0);
//                
//                ExactSol_Laplacian_Ex4(x, sol0, deriv);
//                TPZConnect * con = &(celskeleton->Connect(0));
//                int seq = con->SequenceNumber();
//                cmesh->Block()(seq,0,0,0) = sol0[0];
//                std::cout << sol0[0] << std::endl;
//                
//                x[0] = coord(0,1);
//                x[1] = coord(1,1);
//                ExactSol_Laplacian_Ex4(x, sol1, deriv);
//                con = &(celskeleton->Connect(1));
//                seq = con->SequenceNumber();
//                cmesh->Block()(seq,0,0,0) = sol1[0];
//                std::cout << sol1[0] << std::endl;
//                
//                con = &(celskeleton->Connect(2));
//                seq = con->SequenceNumber();
//                cmesh->Block()(seq,0,0,0) = sol1[0] + sol0[0] - h*h/4.;
//                std::cout << sol1[0] + sol0[0] - h*h/4. << std::endl;
//                
//            }
//            if (coord(1,0) != coord(1,1)) {
//                TPZFNMatrix<200,std::complex<double>> coef(3,1,0);
//                TPZManVector<REAL> x(2,0);
//                TPZManVector<REAL> sol(1,0);
//                TPZFMatrix<STATE> deriv;
//                
//                x[0] = coord(0,0);
//                x[1] = coord(1,0);
//                
//                ExactSol_Laplacian_Ex4(x, sol, deriv);
//                
//                TPZConnect * con = &(celskeleton->Connect(0));
//                int seq = con->SequenceNumber();
//                cmesh->Block()(seq,0,0,0) = sol[0];
//                
//                con = &(celskeleton->Connect(2));
//                seq = con->SequenceNumber();
//                cmesh->Block()(seq,0,0,0) = sol[0];
//                
//                con = &(celskeleton->Connect(1));
//                seq = con->SequenceNumber();
//                cmesh->Block()(seq,0,0,0) = sol[0];
//                std::cout << sol[0] << std::endl;
//            }
        }
        celelem->Coeficients().Print("solel = ", std::cout, EMathematicaInput);
//        celelem->LoadSolution();
//        TPZVec<REAL> errors;
//        
//        TPZCondensedCompEl * condcel = new TPZCondensedCompEl(cel);
//        condcel->LoadElementReference();
//        TPZElementMatrix ek;
//        TPZElementMatrix ef;
//        condcel->CalcStiff(ek, ef);
//        condcel->LoadSolution();
//        
//        celelem->EvaluateError(ExactSol_Laplacian_Ex4, errors, false);
//        std::stringstream sout;
//        sout << "Scalar_ElementError.csv";
//        
//        
//        std::ofstream results(sout.str(),std::ios::app);
//        for (int ier=0; ier<errors.size(); ier++) {
//            results << errors[ier] << ", ";
//        }
//        results << cmesh->GetDefaultOrder() << ", " << h << endl;
        
    }
}
