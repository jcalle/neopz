
#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "Common.h"
#include "TPZSBFemElementGroup.h"
#include <getopt.h>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.sbfem"));
//static LoggerPtr loggerBF(Logger::getLogger("pz.mesh.sbfemelementgroupBF"));
#endif

void IntegrateDirect(TPZCompMesh *cmesh);
void forcefunction(const TPZManVector<REAL> &co, TPZManVector<REAL> &result);
void solexact(const TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix<STATE> &deriv);

/*
 Example 1: u(x,y) = 1./2*(y*y-x*x), Laplaciano = 0;
 Example 2: u(x,y) = 1/2 * Sin(Pi*x)*Sinh(Pi*y), Laplaciano = 0;
 Example 3: u(x,y) = 1/2 * (Sin(Pi*x)*Sinh(Pi*y) - x*x), Laplaciano = 1;
 Example 4: u(x,y) = (1-x*x)/2, Laplaciano = 1;
 Example 5: u(x,y) = Cos(Pi*x/2)*Cos(Pi*y/2)
 
 */

int main(int argc, char *argv[])
{
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
    int option_index = 0;
    struct option sbfem_options[]{
        {"scalarproblem", required_argument, 0, 's'}, // name, has_arg, flag, val
        {"bodyforces", required_argument, 0, 'b'},
        {"maxnelx", required_argument, 0, 'm'},
        {"porder", required_argument, 0, 'p'},
        {"internalporder", required_argument, 0, 'i'},
        {"usesbfem", required_argument, 0,'u'},
        {"usepolynomialfunctions", required_argument, 0, 'y'},
        {"example", required_argument, 0, 'e'},
        {"shapefunctionsplot", no_argument, 0, 'f'},
        {0,0,0,0}
    };
    int c;
//    bool bodyforces = 0;
//    int maxnelxcount = 4;
//    bool plotshapefunctions = 0;
    bool useexact = true;
    
    int numrefskeleton = 0;
    
    bool usesbfem = true, scalarproblem = true, bodyforces = false, plotshapefunctions = false, usepoly = false;
    int maxnelxcount, maxporder, example, maxinternalporder;
    
    while ((c = getopt_long(argc, argv, "s:b:m:p:u:y:f", sbfem_options, &option_index)) != -1){
        switch (c) {
            case 's':
                scalarproblem = atoi(optarg);
                break;
            case 'b':
                bodyforces = atoi(optarg);
                break;
            case 'm':
                maxnelxcount = atoi(optarg);
                break;
            case 'p':
                maxporder = atoi(optarg);
                break;
            case 'u':
                usesbfem = atoi(optarg);
                break;
            case 'i':
                maxinternalporder = atoi(optarg);
                break;
            case 'y':
                usepoly = atoi(optarg);
                break;
            case 'e':
                example = atoi(optarg);
                break;
            case 'f':
                plotshapefunctions = false;
                break;
            case 0:
                break;
            default:
                break;
        }
    }
    if (usesbfem == false) {
        numrefskeleton = 1;
    }
    
    int counter = 1;
    
    usesbfem = true;
    scalarproblem = false;
    bodyforces = true;
    usepoly = false;
    maxporder = 7;
    maxinternalporder = 2;
    maxnelxcount = 4;
    example = 5;
    SetExample(example);
    
#ifdef _AUTODIFF
    LaplaceExact.fExact = TLaplaceExample1::ECosCos;
//    LaplaceExact.fExact = TLaplaceExample1::ESinSin;
//    LaplaceExact.fExact = TLaplaceExample1::E10SinSin;
//    LaplaceExact.fExact = TLaplaceExample1::ESinCos;
//    LaplaceExact.fExact = TLaplaceExample1::EExactTest;
//    ElastExact.fProblemType = TElasticity2DAnalytic::EDispx;
//    ElastExact.fProblemType = TElasticity2DAnalytic::EDispy;
//    ElastExact.fProblemType = TElasticity2DAnalytic::ERot;
//    ElastExact.fProblemType = TElasticity2DAnalytic::EShear;
//    ElastExact.fProblemType = TElasticity2DAnalytic::ELoadedBeam;
//    ElastExact.fProblemType = TElasticity2DAnalytic::EPoly;
    ElastExact.fProblemType = TElasticity2DAnalytic::ESin;
    ElastExact.fE = 2.6;
    ElastExact.fPoisson = 0.3;
//    ElastExact.fPlaneStress = 1;
#endif
    
    std::stringstream sout;
    if (scalarproblem)
    {
        sout << "ScalarRegularSolution.csv";
    }
    else{
        sout << "Elastic2DRegularSolution.csv";
    }
    TPZManVector<std::string> description(12);
    description[0] = "Example";
    description[1] = "p";
    description[2] = "pinternal";
    description[3] = "h";
    description[4] = "Number of equations";
    description[5] = "Error Energy";
    description[6] = "Error L2";
    description[7] = "Error Seminorm Energy";
    description[8] = "Scalar problem?";
    description[9] = "With SBFem?";
    description[10] = "With internal polynomials?";
    description[11] = "With polynomial functions?";
    
    std::ofstream descr(sout.str(),std::ios::app);
    
    for (int64_t i=0; i<description.size()-1; i++) {
        descr << description[i] << ", ";
    }
    descr << description[description.size()-1] << std::endl;
    
    for ( int POrder = 1; POrder <= maxporder; POrder += 1)
    {
        int pinternal = POrder;
//        for (int pinternal = 1; pinternal <= maxinternalporder; pinternal++)
        {
            if (POrder == 3 && !scalarproblem) {
                maxnelxcount = 3;
            }
            for(int nelxcount = 1; nelxcount <= maxnelxcount; nelxcount++)
            {
                int nelx = 1 << (nelxcount-1);
                bool useexact = true;
                
                if (bodyforces) {
                    TPZSBFemElementGroup::gDefaultPolynomialOrder = pinternal;
                }
                if (usepoly) {
                    TPZSBFemElementGroup::gPolynomialShapeFunctions = usepoly;
                }
                
                TPZCompMesh *SBFem;
                if(usesbfem)
                {
                    SBFem = SetupSquareMesh(nelx,numrefskeleton,POrder, scalarproblem,useexact);
                }
                else
                {
                    SBFem = SetupSquareH1Mesh(nelx, POrder, scalarproblem, useexact);
                }
//                SBFem->Print();
                
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
                std::cout << "irefskeleton = " << numrefskeleton << std::endl;
                std::cout << "POrder = " << POrder << std::endl;
                
                bool mustOptimizeBandwidth = false;
                TPZAnalysis * Analysis = new TPZAnalysis(SBFem,mustOptimizeBandwidth);
                Analysis->SetStep(counter++);
                std::cout << "neq = " << SBFem->NEquations() << std::endl;
                
                SolveSist(Analysis, SBFem);
                
                std::cout << "Post processing\n";
//    #ifdef _AUTODIFF
//                if(scalarproblem)
//                {
//                    //Analysis->SetExact(Laplace_exact);
//                    Analysis->SetExact(ExactSol_Laplacian_Ex4);
//                }
//                else
//                {
//                    Analysis->SetExact(Elasticity_exact);
//                }
//    #endif
                if(scalarproblem){
                    switch (example) {
                        case 1:
                            Analysis->SetExact(ExactSol_Laplacian_Ex1);
                            break;
                            
                        case 2:
                            Analysis->SetExact(ExactSol_Laplacian_Ex2);
                            break;
                            
                        case 3:
                            Analysis->SetExact(ExactSol_Laplacian_Ex3);
                            break;
                            
                        case 4:
                            Analysis->SetExact(ExactSol_Laplacian_Ex4);
                            break;
                            
                        case 5:
                            Analysis->SetExact(Laplace_exact);
                            break;
                            
                        default:
                            break;
                    }
                }
                else{
                    Analysis->SetExact(Elasticity_exact);
//                    Analysis->SetExact(ExactSol_ElasticityOoietal1);
                }
                
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
                
                if(1)
                {
                    std::ofstream out("../CompMeshWithSol.txt");
                    SBFem->Print(out);
                }
                
                TPZFNMatrix<3,REAL> sol0 = SBFem->Solution();
                if(plotshapefunctions)
                {
                    for (int i=0; i<sol0.Rows() ;i++){
                        
                        TPZFNMatrix<3,REAL> sol = SBFem->Solution();
                        sol.Zero();
                        sol(i,0) = 1;
                        
                        SBFem->LoadSolution(sol);
                        Analysis->LoadSolution(sol);
                        
                        TPZStack<std::string> vecnames,scalnames;
                        // scalar
                        scalnames.Push("State");
                        Analysis->DefineGraphMesh(2, scalnames, vecnames, "../ShapeFunctions.vtk");
                        Analysis->PostProcess(3);
                    }
                    SBFem->LoadSolution(sol0);
                    Analysis->LoadSolution(sol0);
                }
                
//                sol0.Print("fcoef = ", std::cout, EMathematicaInput);

                std::cout << "Compute errors\n";
                
                Analysis->SetThreadsForError(8);
                TPZVec<REAL> errors(3,0);
                bool store_error = false;
                Analysis->PostProcessError(errors,store_error);
                
                std::ofstream results(sout.str(),std::ios::app);
                results << example << ", " << POrder << ", " << pinternal << ", " << 2./nelx << ", " << SBFem->NEquations() << ", ";
                //            for(int i=0;i<3;i++) errors[i] *= 1e6;
                for (int ier=0; ier<errors.size(); ier++) {
                    results << errors[ier] << ", ";
                }
                if (usesbfem) {
                    results << "Yes" << ", ";
                } else{
                    results << "No" << ", ";
                }
                if (scalarproblem) {
                    results << "Yes" << ", ";
                } else{
                    results << "No" << ", ";
                }
                if (bodyforces) {
                    results << "Yes" << ", ";
                } else{
                    results << "No" << ", ";
                }
                if (usepoly) {
                    results << "Yes" << std::endl;
                } else{
                    results << "No" << std::endl;
                }
                
//                REAL h = 2./REAL(nelx);
//                SetInterpolation(SBFem, h);
//                if(scalarproblem)
//                {
//                    TPZStack<std::string> vecnames,scalnames;
//                    // scalar
//                    scalnames.Push("State");
//                    Analysis->DefineGraphMesh(2, scalnames, vecnames, "../RegularSolution.vtk");
//                    Analysis->PostProcess(3);
//                }
//                Analysis->PostProcessError(errors,store_error);
//                SBFem->Print(std::cout);
                
                delete Analysis;
                delete SBFem;
            }
        }
        
    }
    std::cout << "Check:: Calculation finished successfully" << std::endl;
    return EXIT_SUCCESS;
}
