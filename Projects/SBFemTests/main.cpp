#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "Common.h"
#include "TPZSBFemElementGroup.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.sbfem"));
//static LoggerPtr loggerBF(Logger::getLogger("pz.mesh.sbfemelementgroupBF"));
#endif

void IntegrateDirect(TPZCompMesh *cmesh);
void forcefunction(const TPZManVector<REAL> &co, TPZManVector<REAL> &result);
void solexact(const TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix<STATE> &deriv);

int main(int argc, char *argv[])
{
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    bool scalarproblem = true;
    
    int maxnelxcount = 6;
    int numrefskeleton = 1;
    int maxporder = 2;
    int counter = 1;
    bool usesbfem = true;
    if (usesbfem == false) {
        numrefskeleton = 1;
    }
#ifdef _AUTODIFF
        ElastExact.fProblemType = TElasticity2DAnalytic::EBend;
        LaplaceExact.fExact = TLaplaceExample1::ESinSin;
        LaplaceExact.fExact = TLaplaceExample1::E10SinSin; //NÃO FUNCIONOU P/ P=2
        LaplaceExact.fExact = TLaplaceExample1::ECosCos;
        LaplaceExact.fExact = TLaplaceExample1::ESinCos;
        LaplaceExact.fExact = TLaplaceExample1::EExactTest;
#endif
    for ( int POrder = 2; POrder <= maxporder; POrder += 1)
    {
        for (int irefskeleton = 0; irefskeleton < numrefskeleton; irefskeleton++)
        {
            if (POrder == 3 && !scalarproblem) {
                maxnelxcount = 3;
            }
            for(int nelxcount = 1; nelxcount <= maxnelxcount; nelxcount++)
            {
                int nelx = 1 << (nelxcount-1);
                bool useexact = true;
                if(!scalarproblem)
                {
#ifdef _AUTODIFF
                    ElastExact.fE = 10;
                    ElastExact.fPoisson = 0.3;
                    ElastExact.fPlaneStress = 0;
#endif
                }
                TPZSBFemElementGroup::gDefaultPolynomialOrder = POrder;
                TPZCompMesh *SBFem;
                if(usesbfem)
                {
                    SBFem = SetupSquareMesh(nelx,irefskeleton,POrder, scalarproblem,useexact);
                }
                else
                {
                    SBFem = SetupSquareH1Mesh(nelx, POrder, scalarproblem, useexact);
                }
                if(0 && !scalarproblem)
                {
#ifdef _AUTODIFF
                    ElastExact.fProblemType = TElasticity2DAnalytic::EBend;
                    TPZManVector<REAL,3> x(3,0.);
                    TPZFNMatrix<4,STATE> tensor(2,2);
                    for(int i=-1; i<3; i+=2)
                    {
                        for (int j=-1; j<3; j+=2) {
                            x[0] = i;
                            x[1] = j;
                            ElastExact.Sigma(x, tensor);
                            std::cout << "x = " << x << " tensor " << tensor << std::endl;
                        }
                    }
#endif
                }
                SBFem->ComputeNodElCon();
                SBFem->InitializeBlock();
#ifdef LOG4CXX
                if(logger->isDebugEnabled())
                {
                    std::stringstream sout;
                    SBFem->Print(sout);
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
                
                std::cout << "nelx = " << nelx << std::endl;
                std::cout << "irefskeleton = " << irefskeleton << std::endl;
                std::cout << "POrder = " << POrder << std::endl;
                
                // Visualization of computational meshes
                
                bool mustOptimizeBandwidth = true;
                TPZAnalysis * Analysis = new TPZAnalysis(SBFem,mustOptimizeBandwidth);
                Analysis->SetStep(counter++);
                std::cout << "neq = " << SBFem->NEquations() << std::endl;
                SolveSist(Analysis, SBFem);
                
                std::cout << "Post processing\n";
                //        ElasticAnalysis->Solution().Print("Solution");
                //        mphysics->Solution().Print("expandec");
                //                Analysis->SetExact(ExactSol_Laplacian);
#ifdef _AUTODIFF
                if(scalarproblem)
                {
                    Analysis->SetExact(ExactSol_Laplacian);
                    //                    Analysis->SetExact(solexact);
                    //                    Analysis->SetExact(Laplace_exact);
                }
                else
                {
                    Analysis->SetExact(Elasticity_exact);
                }
#endif
                //                ElasticAnalysis->SetExact(Singular_exact);
                //                std::cout << SBFem->Solution() << std::endl;
                int64_t neq = SBFem->Solution().Rows();
                
                if(scalarproblem)
                {
                    TPZStack<std::string> vecnames,scalnames;
                    // scalar
                    scalnames.Push("State");
                    Analysis->DefineGraphMesh(2, scalnames, vecnames, "../RegularSolution.vtk");
                    Analysis->PostProcess(3);
                }
                else
                {
                    TPZStack<std::string> vecnames,scalnames;
                    // scalar
                    vecnames.Push("Displacement");
                    scalnames.Push("SigmaX");
                    scalnames.Push("SigmaY");
                    scalnames.Push("TauXY");
                    scalnames.Push("EpsX");
                    scalnames.Push("EpsY");
                    scalnames.Push("EpsXY");
                    std::stringstream sout;
                    sout << "../RegularElasticity2DSolution";
                    if(usesbfem)
                    {
                        sout << "_SBFem.vtk";
                    }
                    else
                    {
                        sout << "_H1.vtk";
                    }
                    Analysis->DefineGraphMesh(2, scalnames, vecnames,sout.str());
                    Analysis->PostProcess(3);
                }
                
                if(0)
                {
                    std::ofstream out("../CompMeshWithSol.txt");
                    SBFem->Print(out);
                }
                
                std::cout << "Compute errors\n";
                
                TPZManVector<REAL,10> errors(3,0.);
                Analysis->SetThreadsForError(8);
                TPZVec<REAL> error(3,0);
                Analysis->PostProcessError(error);
                
                delete Analysis;
                delete SBFem;
            }
        }
    }
    std::cout << "Check:: Calculation finished successfully" << std::endl;
    return EXIT_SUCCESS;
}

void UniformRefinement(TPZGeoMesh *gMesh, int nh)
{
    for ( int ref = 0; ref < nh; ref++ ){
        TPZVec<TPZGeoEl *> filhos;
        int64_t n = gMesh->NElements();
        for ( int64_t i = 0; i < n; i++ ){
            TPZGeoEl * gel = gMesh->ElementVec() [i];
            if (gel->Dimension() == 2 || gel->Dimension() == 1) gel->Divide (filhos);
        }//for i
    }//ref
}

#include "TPZSBFemVolume.h"
#include "TPZSBFemElementGroup.h"

void IntegrateDirect(TPZCompMesh *cmesh)
{
    int64_t nel = cmesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZSBFemElementGroup *elgr = dynamic_cast<TPZSBFemElementGroup *>(cel);
        if (elgr) {
            TPZStack<TPZCompEl *,5> elstack = elgr->GetElGroup();
            int nvol = elstack.size();
            TPZElementMatrix ekvol, efvol, ekgrp, efgrp;
            elgr->CalcStiff(ekgrp, efgrp);
            for (int iv=0; iv<nvol; iv++) {
                TPZCompEl *vcel = elstack[iv];
                TPZSBFemVolume *elvol = dynamic_cast<TPZSBFemVolume *>(vcel);
                TPZElementMatrix ek,ef;
                elvol->CalcStiff(ek, ef);
                if (iv==0) {
                    ekvol = ek;
                    efvol = ef;
                }
                else
                {
                    ekvol.fMat += ek.fMat;
                    efvol.fMat += ef.fMat;
                }
            }
            //            ekgrp.fMat.Print("EKGRP = ",std::cout,EMathematicaInput);
            //            ekvol.fMat.Print("EKVOL = ",std::cout,EMathematicaInput);
            ekvol.fMat -= ekgrp.fMat;
            std::cout << "IntegrateDirect Norm of difference " << Norm(ekvol.fMat) << std::endl;
            break;
        }
    }
    
    
}
