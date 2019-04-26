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

#ifdef _AUTODIFF
TElasticity2DAnalytic ElastExact;

TLaplaceExample1 LaplaceExact;

#endif

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.sbfem"));
//static LoggerPtr loggerBF(Logger::getLogger("pz.mesh.sbfemelementgroupBF"));
#endif

//TElasticity2DAnalytic::EDefState TElasticity2DAnalytic::fProblemType = TElasticity2DAnalytic::EStretchx;

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

void Constant(const TPZVec<REAL> &x, TPZVec<STATE> &val)
{
    val.resize(1);
    val[0] = 1.;
}

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
//    int dim = 2;

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
        REAL lamelambda = 0.,lamemu = 0.5e3, fx= 0, fy = 0.5;
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
        
        int polynomialorder = 2;
//        matloc->SetForcingFunction(LaplaceExact.ForcingFunction());
        TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(Constant, polynomialorder);
        TPZAutoPointer<TPZFunction<STATE> > autodummy = dummy;
        matloc->SetForcingFunction(autodummy);
    
//        TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(forcefunction, polynomialorder);
//        TPZAutoPointer<TPZFunction<STATE> > autodummy = dummy;
//        matloc->SetForcingFunction(autodummy);
    
        material = matloc;
        nstate = 1;
    }
    
    
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,2,0.);    val2(0,0) = 0.0;
    
    TPZMaterial * BSkeleton = material->CreateBC(material,ESkeleton,1, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(material,Ebc2, 0, val1, val2);
    TPZMaterial * BCond4 = material->CreateBC(material,Ebc4, 0, val1, val2);
    
    val1(0,0) = 1.; val2(0,0) = 1.;
    TPZMaterial * BCond1 = material->CreateBC(material,Ebc1, 1, val1, val2);
//    TPZMaterial * BCond2 = material->CreateBC(material,Ebc2, 1, val1, val2);
    TPZMaterial * BCond3 = material->CreateBC(material,Ebc3, 1, val1, val2);
//    TPZMaterial * BCond4 = material->CreateBC(material,Ebc4, 1, val1, val2);
//    TPZMaterial * BCond1 = material->CreateBC(material,Ebc1, 0, val1, val2);
//    TPZMaterial * BCond2 = material->CreateBC(material,Ebc2, 0, val1, val2);
//    TPZMaterial * BCond3 = material->CreateBC(material,Ebc3, 0, val1, val2);
//    TPZMaterial * BCond4 = material->CreateBC(material,Ebc4, 0, val1, val2);
    
    
    int polynomialorder = 2;
    //        matloc->SetForcingFunction(LaplaceExact.ForcingFunction());
    TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(Constant, polynomialorder);
    TPZAutoPointer<TPZFunction<STATE> > autodummy = dummy;
//    BCond1->SetBCForcingFunction(LaplaceExact.TensorFunction());
    BCond2->SetBCForcingFunction(autodummy);
//    BCond3->SetBCForcingFunction(LaplaceExact.TensorFunction());
    BCond4->SetBCForcingFunction(autodummy);
    
    cmesh->InsertMaterialObject(material);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
    cmesh->InsertMaterialObject(BSkeleton);
    
    val1.Zero();
    val2.Zero();
    TPZMaterial * BCond5 = material->CreateBC(material,EBCPoint1, 0, val1, val2);
    TPZMaterial * BCond6 = material->CreateBC(material,EBCPoint2, 0, val1, val2);
    cmesh->InsertMaterialObject(BCond5);
    cmesh->InsertMaterialObject(BCond6);
    
}

TPZCompMesh *SetupSquareMesh(int nelx, int nrefskeleton, int porder, bool scalarproblem, bool useexact)
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
    gengrid.Read(gmesh, EGroup);
    gengrid.SetBC(gmesh, 4, Ebc1);
    gengrid.SetBC(gmesh, 5, Ebc2);
    gengrid.SetBC(gmesh, 6, Ebc3);
    gengrid.SetBC(gmesh, 7, Ebc4);
//    {
//        TPZManVector<int64_t,2> nodeindex(1);
//		int64_t index;
//        nodeindex[0] = 0;
//        gmesh->CreateGeoElement(EPoint, nodeindex, EBCPoint1, index);
//        nodeindex[0] = nelx;
//        gmesh->CreateGeoElement(EPoint, nodeindex, EBCPoint2, index);
//        gmesh->BuildConnectivity();
//    }
    
    std::map<int,int> matmap;
    matmap[EGroup] = 1;
    TPZBuildSBFem build(gmesh,ESkeleton,matmap);
    
    build.StandardConfiguration();
    
//    nrefskeleton=0;
    build.DivideSkeleton(nrefskeleton);
//    AddSkeletonElements(gmesh);
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
    InsertMaterialObjects(SBFem,!elasticityproblem, useexact);
    
    int polynomialorder=1;
//    if(problemtype == 1)
//    {
//        TPZMaterial *BCond2 = SBFem->FindMaterial(Ebc2);
//        TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(HarmonicNeumannRight, polynomialorder);
//        TPZAutoPointer<TPZFunction<STATE> > autodummy = dummy;
//        BCond2->SetForcingFunction(autodummy);
//    }
//    if(problemtype == 1)
//    {
//        TPZMaterial *BCond4 = SBFem->FindMaterial(Ebc4);
//        TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(HarmonicNeumannLeft, polynomialorder);
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
    
    
    {
        std::ofstream gout("Geometry_SBFEM.vtk");
        TPZVTKGeoMesh vtk;
        vtk.PrintGMeshVTK(gmesh, gout, true);
        std::cout << Ebc1 << ", " << Ebc2 << ", " << Ebc3 << ", " << Ebc4 << std::endl;
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
//    {
//        TPZManVector<int64_t,2> nodeindex(1);
//		int64_t index;
//        nodeindex[0] = 0;
//        gmesh->CreateGeoElement(EPoint, nodeindex, EBCPoint1, index);
//        nodeindex[0] = nelx;
//        gmesh->CreateGeoElement(EPoint, nodeindex, EBCPoint2, index);
//    }
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
    InsertMaterialObjects(SBFem,!elasticityproblem, useexact);
//    int polynomialorder=1;
//    if(problemtype == 1)
//    {
//        TPZMaterial *BCond2 = SBFem->FindMaterial(Ebc2);
//        TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(HarmonicNeumannRight,polynomialorder);
//        TPZAutoPointer<TPZFunction<STATE> > autodummy = dummy;
//        BCond2->SetForcingFunction(autodummy);
//    }
//    if(problemtype == 1)
//    {
//        TPZMaterial *BCond4 = SBFem->FindMaterial(Ebc4);
//        TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(HarmonicNeumannLeft, polynomialorder);
//        TPZAutoPointer<TPZFunction<STATE> > autodummy = dummy;
//        BCond4->SetForcingFunction(autodummy);
//    }
    

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


void Constant_Laplacian(const TPZVec<REAL> &x, TPZVec<STATE> &val)
{
    val[0] = 1.0;
    val[1] = 0;
}

void ExactSol_Laplacian(const TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix<STATE> &deriv){
//    REAL Pi = M_PI;
//    REAL xi = x[0];
//    REAL eta = x[1];
//    sol.resize(1);
//    
//    sol[0] = (32*(6043467208275*cos((Pi*eta)/2.)*(839475*cos((Pi*xi)/2.) + 12915*cos((5*Pi*xi)/2.) +
//                                            2275*cos((9*Pi*xi)/2.) + 55965*sin((3*Pi*(1 + xi))/2.) +
//                                            4797*sin((7*Pi*(1 + xi))/2.)) +
//         297305883*cos((5*Pi*eta)/2.)*(262528875*cos((Pi*xi)/2.) +
//                                          27303003*cos((5*Pi*xi)/2.) +
//                                          1625*(4403*cos((9*Pi*xi)/2.) + 41181*sin((3*Pi*(1 + xi))/2.) +
//                                                8109*sin((7*Pi*(1 + xi))/2.))) +
//         25*(58671*(6588266139*cos((Pi*xi)/2.)*
//                    (35*sin((3*Pi*(1 + eta))/2.) + 3*sin((7*Pi*(1 + eta))/2.)) +
//                    52370955*cos((5*Pi*xi)/2.)*
//                    (259*sin((3*Pi*(1 + eta))/2.) + 51*sin((7*Pi*(1 + eta))/2.)) +
//                    3145*(9947*sin((9*Pi*xi)/2.)*
//                          (91*sin((3*Pi*(1 + eta))/2.) + 27*sin((7*Pi*(1 + eta))/2.)) +
//                          66885*sin((3*Pi*(1 + eta))/2.)*
//                          (203*sin((3*Pi*(1 + xi))/2.) + 27*sin((7*Pi*(1 + xi))/2.)) +
//                          5265*sin((7*Pi*(1 + eta))/2.)*
//                          (343*sin((3*Pi*(1 + xi))/2.) + 87*sin((7*Pi*(1 + xi))/2.)))) +
//             31283315*cos((9*Pi*eta)/2.)*
//             (17579835*cos((Pi*xi)/2.) +
//              41*(66339*cos((5*Pi*xi)/2.) +
//                  53*(455*cos((9*Pi*xi)/2.) + 2457*sin((3*Pi*(1 + xi))/2.) +
//                      729*sin((7*Pi*(1 + xi))/2.)))))))/(5.073339634666656e18*pow(Pi,4));
    sol[0] = 1/2*(x[0]*x[0]-1);
}

void forcefunction(const TPZVec<REAL> &co, TPZVec<STATE> &result)
{
    result.resize(1);
    
    double xv = co[0];
    double yv = co[1];
    
    double Pi = M_PI;
    result[0] = 8.0*Pi*Pi*cos(2.0*Pi*yv)*sin(2.0*Pi*xv);
//    result[1] = -8.0*Pi*Pi*cos(2.0*Pi*yv)*sin(2.0*Pi*xv);
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
