/**
 * @file
 * @brief Performs modal analysis in an electromagnetic waveguide.
 * @details This project is aimed at electromagnetic analysis
 * of waveguides, but it is not restricted to it. More generally,
 * it could be a starting point for problems with longitudinal axis symmetry.
 * The problem is stated as a generalised eigenvalue problem, as seen in
 * Jin, Jian-Ming. The finite element method in electromagnetics.(contd.)
 * John Wiley & Sons, 2015 (Chapter 8).
 * The transverse field components are approximated using H(curl,\Omega)
 * functions,
 * and the longitudinal field components are approximated using H^1(\Omega)
 * functions.
 *
 * @author Francisco Orlandini
 * @since 2015
 */

#include <stddef.h>               // for NULL
#include <fstream>
#include <TPZEigenAnalysis.h>
#include <TPZSpStructMatrix.h>
#include <tpzgeoelmapped.h>
#include <tpzarc3d.h>
#include "pzgengrid.h"
#include "pzgmesh.h"
#include "TPZVTKGeoMesh.h"
#include "pzbndcond.h"
#include "pzlog.h"
#include "pzstrmatrix.h"
#include "pzl2projection.h"
#include "pzfstrmatrix.h"
#include "pzbuildmultiphysicsmesh.h"
#include "boost/date_time/posix_time/posix_time.hpp"
#ifdef USING_SLEPC
#include <TPZSlepcEPSHandler.h>
#include <TPZSlepcSTHandler.h>
#include <Complex/TPZMatWaveguidePml.h>
#include <tpzgeoblend.h>

#endif
#include "parameter_handler.h"
#include "TPZMatWaveguideCutOffAnalysis.h"
#include "TPZMatModalAnalysis.h"
#include "pzintel.h"
#include "TPZGmshReader.h"
#include "SPZModalAnalysisDataReader.h"
#include "SPZModalAnalysisData.h"

void RunSimulation(SPZModalAnalysisData &simData,std::ostringstream &eigeninfo, const int &factor);

typedef struct QuadrilateralData QuadrilateralData;

void
CreateGMeshRectangularWaveguide(TPZGeoMesh *&gmesh, const std::string mshFileName, TPZVec<int> &matIdVec,
                                const std::string &prefix, const bool &print, const REAL &scale, const int &factor,
                                const REAL wDomain, const REAL hDomain, TPZVec<SPZModalAnalysisData::boundtype> &boundTypeVec,
                                TPZVec<SPZModalAnalysisData::pmltype> &pmlTypeVec, const bool &usingSymmetry,
                                const SPZModalAnalysisData::boundtype &symmetryType);
void
CreateStepFiberMesh(TPZGeoMesh *&gmesh, const std::string mshFileName, TPZVec<int> &matIdVec, const std::string &prefix,
                    const bool &print, const REAL &scale, const int &factor,
                    TPZVec<SPZModalAnalysisData::boundtype> &boundTypeVec,
                    TPZVec<SPZModalAnalysisData::pmltype> &pmlTypeVec, const REAL &realRCore, const REAL &dPML,
                    const REAL &boundDist, const REAL &outerReffIndex);

void
ReadGMesh(TPZGeoMesh *&gmesh, const std::string mshFileName, TPZVec<int> &matIdVec, const std::string &prefix, const bool &print, const REAL &scale = 1.);

void CreateCMesh(TPZVec<TPZCompMesh *> &meshVecOut, TPZGeoMesh *gmesh, int pOrder, TPZVec<int> &matIdVec,
                 TPZVec<STATE> &urVec, TPZVec<STATE> &erVec, REAL f0, bool isCutOff, const std::string &prefix,
                 const bool &print, const REAL &scale, const REAL &alphaMax,
                 TPZVec<SPZModalAnalysisData::boundtype> &boundTypeVec,
                 TPZVec<SPZModalAnalysisData::pmltype> &pmlTypeVec);

void FilterBoundaryEquations(TPZVec<TPZCompMesh *> cmeshMF,
                             TPZVec<long> &activeEquations, int &neq,
                             int &neqOriginal);


void CreateGmshMesh(const std::string &meshName, const std::string &newName, const REAL &factor, const int &nThreads,
                    const REAL &scale, const int &meshOrder, const REAL &refIndexOuterDomain);

int main(int argc, char *argv[]) {
#ifdef LOG4CXX
    InitializePZLOG();
#endif



    ParameterHandler prm;
    SPZModalAnalysisDataReader reader(prm,argc,argv);
    SPZModalAnalysisData simData;
    reader.ReadParameters(simData);

    const std::string meshOriginal = simData.pzOpts.meshFile;
    const int pOrderOrig = simData.pzOpts.pOrder;

    boost::posix_time::ptime t1_total =
            boost::posix_time::microsec_clock::local_time();
    std::string eigenFileName = simData.pzOpts.prefix;
    if(simData.pzOpts.exportEigen){
        eigenFileName +="mapord_"+std::to_string(simData.pzOpts.meshOrder)+"_";
        if(simData.physicalOpts.freqVec.size() > 1){
            eigenFileName+="from_l_"+std::to_string(simData.physicalOpts.freqVec[0]);
            eigenFileName+="_to_"+std::to_string(simData.physicalOpts.freqVec[simData.physicalOpts.freqVec.size()])+"_";
        }else{
            eigenFileName+="l_"+std::to_string(simData.physicalOpts.freqVec[0])+"_";
        }
        if(simData.pzOpts.pSteps>1){
            eigenFileName+="from_p_"+std::to_string(pOrderOrig)+"_to_";
            eigenFileName+=std::to_string(pOrderOrig+simData.pzOpts.pSteps-1)+"_";
        }else{
            eigenFileName+="p_"+std::to_string(pOrderOrig);
        }
        eigenFileName+="_n_hsteps_"+std::to_string(simData.pzOpts.hSteps)+".csv";

        std::ofstream file(eigenFileName.c_str(), std::ios::trunc);
        file.close();
    }
    for (int iFreq = 0; iFreq < simData.physicalOpts.freqVec.size(); ++iFreq) {
        simData.physicalOpts.lambda = simData.physicalOpts.freqVec[iFreq];
        simData.physicalOpts.lambda = simData.physicalOpts.isLambda ? simData.physicalOpts.lambda : 299792458 / simData.physicalOpts.lambda;
        simData.pzOpts.scaleFactor = simData.pzOpts.scaleByk0 ? 2 * M_PI / simData.physicalOpts.lambda : simData.pzOpts.scaleFactor;
        if(simData.physicalOpts.freqVec.size() > 1){
            std::cout<<"Beginning step "<<iFreq+1<<" out of "<<simData.physicalOpts.freqVec.size()<<" freq steps."<<std::endl;
            std::cout<<"lambda = "<<simData.physicalOpts.lambda<<std::endl;
        }
        for (int iH = 0; iH < simData.pzOpts.hSteps; ++iH) {
            std::cout << "Beginning step " << iH + 1 << " out of " << simData.pzOpts.hSteps << " h steps."
                      << std::endl;
            if(simData.pzOpts.usingNeoPzMesh == false){

                simData.pzOpts.meshFile = meshOriginal.substr(0, meshOriginal.size() - 4);
                if(simData.pzOpts.externGenMesh){
                    if(iH>0){//refine by splitting
                        simData.pzOpts.meshFile += "h" + std::to_string(iH);
                        const std::string lastMesh = iH > 1 ?
                                                     meshOriginal.substr(0, meshOriginal.size() - 4) + "h" + std::to_string(iH-1) + ".msh":
                                                     meshOriginal;
                        const std::string command = "gmsh -v 3 -refine " + lastMesh + " -o " + simData.pzOpts.meshFile + ".msh";
                        std::cout<<"Generating mesh with: "<<std::endl<<command<<std::endl;
                        std::array<char, 128> buffer;
                        std::string result;
                        std::shared_ptr<FILE> pipe(popen(command.c_str(), "r"), pclose);
                        if (!pipe) throw std::runtime_error("popen() failed!");
                        while (fgets(buffer.data(), 128, pipe.get()) != nullptr){
                            result += buffer.data();
                        }
                        std::cout<<result<<std::endl;
                    }
                    simData.pzOpts.meshFile +=".msh";
                }
                else{
                    const REAL factorVal = simData.pzOpts.factorVec[iH];
                    if(simData.physicalOpts.freqVec.size() > 1 ) {
                        simData.pzOpts.meshFile += "f" + std::to_string(iFreq);
                    }
                    simData.pzOpts.meshFile += "ord" + std::to_string(simData.pzOpts.meshOrder);
                    simData.pzOpts.meshFile += "h" + std::to_string(iH);
                    simData.pzOpts.meshFile += ".msh";
                    const int nMaterials = simData.physicalOpts.erVec.size();
                    const REAL refIndexOuterMat =
                            std::sqrt(std::real(simData.physicalOpts.erVec[nMaterials-1]));
                    CreateGmshMesh(meshOriginal, simData.pzOpts.meshFile, factorVal, simData.pzOpts.nThreads,
                                   simData.pzOpts.scaleFactor, simData.pzOpts.meshOrder,
                                   refIndexOuterMat);
                }
            }

            simData.pzOpts.pOrder = pOrderOrig;
            for (int iP = 0; iP < simData.pzOpts.pSteps; ++iP) {
                std::cout<<"freq step: "<<iFreq+1<<" out of "<<simData.physicalOpts.freqVec.size()<<std::endl;
                std::cout<<"h    step: "<< iH+1<<" out of "<<simData.pzOpts.hSteps<<std::endl;
                std::cout<<"Beginning p step "<<iP+1<<" out of "<<simData.pzOpts.pSteps<<std::endl;
                std::ostringstream eigeninfo;
                RunSimulation(simData,eigeninfo, simData.pzOpts.factorVec[iH]);
                if(simData.pzOpts.exportEigen){
                    std::ofstream file(eigenFileName.c_str(),std::ios::app);
                    file << eigeninfo.str();
                    file.close();
                }
                simData.pzOpts.pOrder++;
            }
        }
    }
    boost::posix_time::ptime t2_total =
            boost::posix_time::microsec_clock::local_time();
    std::cout<<"Total time: "<<t2_total-t1_total<<std::endl;
    return 0;
}

void CreateGmshMesh(const std::string &meshName, const std::string &newName, const REAL &factor, const int &nThreads,
                    const REAL &scale, const int &meshOrder, const REAL &refIndexOuterDomain) {
    std::ostringstream str_factor;
    str_factor << std::setprecision(20) << factor;
    std::ostringstream str_scale;
    str_scale<< std::setprecision(20) << scale;
    std::ostringstream str_n_mat;
    str_n_mat<< std::setprecision(20) << refIndexOuterDomain;

    std::string command = "gmsh " + meshName + " -2 -match ";
    command += " -nt " + std::to_string(nThreads);
    command += " -tol 1e-20 ";
    command += " -v 3 ";
    command += " -setnumber scale "+str_scale.str();
    command += " -setnumber factor "+str_factor.str();
    command += " -setnumber nOuterDomain "+str_n_mat.str();
    command += " -order " + std::to_string(meshOrder);
    if( meshOrder > 1 ) command += " -optimize_ho";
    command += " -o " + newName;
    std::cout<<"Generating mesh with: "<<std::endl<<command<<std::endl;

    std::array<char, 128> buffer;
    std::string result;
    std::shared_ptr<FILE> pipe(popen(command.c_str(), "r"), pclose);
    if (!pipe) throw std::runtime_error("popen() failed!");
    while (fgets(buffer.data(), 128, pipe.get()) != nullptr){
        result += buffer.data();
    }
    std::cout<<result<<std::endl;
}

void RunSimulation(SPZModalAnalysisData &simData,std::ostringstream &eigeninfo, const int &factor) {
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    std::cout<<"Creating GMesh..."<<std::endl;
    boost::posix_time::ptime t1_g =
        boost::posix_time::microsec_clock::local_time();
    TPZVec<int> matIdVec;
    TPZVec<SPZModalAnalysisData::boundtype> boundTypeVec(0);
    TPZVec<SPZModalAnalysisData::pmltype> pmlTypeVec(0);
    if(simData.pzOpts.usingNeoPzMesh){
        switch(simData.pzOpts.pzCase){
            case SPZModalAnalysisData::StepFiber:{
                const REAL realRCore = simData.physicalOpts.stepFiberOpts.realRCore;
                const REAL dPML = simData.physicalOpts.stepFiberOpts.dPML;
                const REAL boundDist = simData.physicalOpts.stepFiberOpts.boundDist;
                const REAL outerMaterialReffIndex = std::sqrt(std::real(simData.physicalOpts.erVec[1]));
                CreateStepFiberMesh(gmesh, simData.pzOpts.meshFile, matIdVec, simData.pzOpts.prefix,
                                    simData.pzOpts.exportGMesh, simData.pzOpts.scaleFactor,
                                    factor, boundTypeVec, pmlTypeVec, realRCore, dPML,
                                    boundDist, outerMaterialReffIndex);
                break;
            }
            case SPZModalAnalysisData::RectangularWG:{
                const REAL wDomain = simData.physicalOpts.rectangularWgOpts.wDomain;
                const REAL hDomain = simData.physicalOpts.rectangularWgOpts.hDomain;
                const bool usingSymmetry = simData.physicalOpts.rectangularWgOpts.usingSymmetry;
                const SPZModalAnalysisData::boundtype boundType= simData.physicalOpts.rectangularWgOpts.symmetryType;
                CreateGMeshRectangularWaveguide(gmesh, simData.pzOpts.meshFile, matIdVec, simData.pzOpts.prefix,
                                                simData.pzOpts.exportGMesh,
                                                simData.pzOpts.scaleFactor, factor, wDomain, hDomain, boundTypeVec,
                                                pmlTypeVec,
                                                usingSymmetry, boundType);
                break;
            }
            default:
                DebugStop();
        }
    }else{
        boundTypeVec = simData.physicalOpts.boundTypeVec;
        pmlTypeVec = simData.physicalOpts.pmlTypeVec;
        ReadGMesh(gmesh, simData.pzOpts.meshFile, matIdVec, simData.pzOpts.prefix,
                  simData.pzOpts.exportGMesh);
    }
    boost::posix_time::ptime t2_g =
        boost::posix_time::microsec_clock::local_time();
    std::cout<<"Created!  "<<t2_g-t1_g<<std::endl;
    TPZVec<TPZCompMesh *> meshVec(1);
    std::cout<<"Creating CMesh..."<<std::endl;
    boost::posix_time::ptime t1_c =
        boost::posix_time::microsec_clock::local_time();

    CreateCMesh(meshVec, gmesh, simData.pzOpts.pOrder, matIdVec, simData.physicalOpts.urVec, simData.physicalOpts.erVec,
                simData.physicalOpts.lambda, simData.physicalOpts.isCutOff, simData.pzOpts.prefix,
                simData.pzOpts.exportCMesh, simData.pzOpts.scaleFactor, simData.physicalOpts.alphaMax, boundTypeVec,
                pmlTypeVec); // funcao para criar a malha computacional
    boost::posix_time::ptime t2_c =
        boost::posix_time::microsec_clock::local_time();
    std::cout<<"Created! "<<t2_c-t1_c<<std::endl;
    TPZCompMesh *cmesh = meshVec[0];

    TPZEigenAnalysis an(cmesh, true);

    TPZManVector<long, 1000> activeEquations;
    int neq = 0;
    int neqOriginal = 0;

    TPZAutoPointer<TPZStructMatrix> strmtrx;
    strmtrx = new TPZSpStructMatrix(cmesh);
    strmtrx->SetNumThreads(simData.pzOpts.nThreads);
    FilterBoundaryEquations(meshVec, activeEquations, neq, neqOriginal);
    strmtrx->EquationFilter().SetActiveEquations(activeEquations);
    an.SetStructuralMatrix(strmtrx);

    //const int nSolutions = neq >= simData.solverOpts.eps_nev ? simData.solverOpts.eps_nev : neq;
    TPZSlepcEPSHandler<STATE> solver;
    TPZSlepcSTHandler stHandler;
    solver.SetTolerances(simData.solverOpts.eps_tol,simData.solverOpts.eps_max_its);
    solver.SetConvergenceTest(simData.solverOpts.eps_conv_test);
    solver.SetWhichEigenpairs(simData.solverOpts.eps_which_eig);
    solver.SetTargetEigenvalue(simData.solverOpts.target);

    stHandler.SetPrecond(simData.solverOpts.st_precond);
    stHandler.SetSolver(simData.solverOpts.st_solver);
    stHandler.SetSolverTol(simData.solverOpts.ksp_rtol,simData.solverOpts.ksp_atol,
                           simData.solverOpts.ksp_dtol,simData.solverOpts.ksp_max_its);
    stHandler.SetType(simData.solverOpts.st_type,simData.solverOpts.target);
    solver.SetST(stHandler);

    solver.SetTrueResidual(simData.solverOpts.eps_true_res);
    solver.SetProblemType(simData.solverOpts.eps_prob_type);
    solver.SetType(simData.solverOpts.eps_type);
    solver.SetKrylovOptions(simData.solverOpts.eps_krylov_locking,simData.solverOpts.eps_krylov_restart);
    solver.SetEPSDimensions(simData.solverOpts.eps_nev, simData.solverOpts.eps_ncv, simData.solverOpts.eps_mpd);
    solver.SetVerbose(simData.solverOpts.eps_verbose);
    solver.SetAbsoluteValue(simData.pzOpts.absVal);
    an.SetSolver(solver);

    if(simData.pzOpts.exportGMesh){
        std::cout<<"Printing mesh..."<<std::endl;
        TPZStack<std::string> scalnames, vecnames;
        vecnames.Push("Material");
        std::string plotfile = simData.pzOpts.prefix + "Material" + ".vtk";
        // estara na pasta debug
        const int dim = 2;
        an.DefineGraphMesh(dim, scalnames, vecnames,
                           plotfile);  // define malha grafica
        an.PostProcess(simData.pzOpts.vtkRes);
    }
    std::cout << "Assembling..." << std::endl;
    boost::posix_time::ptime t1 =
        boost::posix_time::microsec_clock::local_time();
    an.Assemble();
    boost::posix_time::ptime t2 =
        boost::posix_time::microsec_clock::local_time();
    std::cout << "Finished assembly." << std::endl;

    std::cout << "Solving..." << std::endl;
    boost::posix_time::ptime t3 =
        boost::posix_time::microsec_clock::local_time();

    an.Solve();
    boost::posix_time::ptime t4 =
        boost::posix_time::microsec_clock::local_time();
    std::cout << "Time for assembly " << t2 - t1 << " Time for solving "
              << t4 - t3 << std::endl;
    const TPZManVector<SPZAlwaysComplex<STATE>::type> eigenValues = an.GetEigenvalues();
    const TPZFMatrix<SPZAlwaysComplex<STATE>::type> eigenVectors = an.GetEigenvectors();

    typedef std::numeric_limits< double > dbl;
    std::cout.precision(dbl::max_digits10);
    for(int i = 0; i < eigenValues.size() ; i++){
        std::cout<<std::fixed<<eigenValues[i]<<std::endl;
    }
    if(simData.pzOpts.exportEigen){
        eigeninfo.precision(dbl::max_digits10);
        std::cout<<"Exporting eigen info..."<<std::endl;
        REAL hSize = -1e12;
        REAL elRadius = 0;
        TPZVec<REAL> qsi(2,0.25);
        TPZVec<REAL> x(3,0.);
        long nel = 0;
        for (int j = 0; j < cmesh->NElements(); ++j) {
            TPZGeoEl &el = *(cmesh->Element(j)->Reference());
            if(el.Type() == MElementType::EQuadrilateral){
                continue;
            }
            nel++;
            TPZBndCond *matBound = dynamic_cast<TPZBndCond *>(meshVec[0]->MaterialVec()[el.MaterialId()]);
            TPZMatWaveguidePml *matPml = dynamic_cast<TPZMatWaveguidePml *>(meshVec[0]->MaterialVec()[el.MaterialId()]);
            if(!matBound && !matPml){
               elRadius = el.ElementRadius();
               hSize = elRadius > hSize ? elRadius : hSize;
            }
        }
        eigeninfo << std::fixed << neq << "," << nel << ",";
        eigeninfo << std::fixed << hSize << "," << simData.pzOpts.pOrder<<",";
        eigeninfo << std::fixed << simData.physicalOpts.lambda<<",";
        eigeninfo << eigenValues.size()<<",";
        for(int i = 0; i < eigenValues.size() ; i++){
            eigeninfo<<std::fixed<<std::real(eigenValues[i])<<",";
            eigeninfo<<std::fixed<<std::imag(eigenValues[i]);
            if(i != eigenValues.size() - 1 ) {
                eigeninfo << ",";
            }
        }
        eigeninfo << std::endl;

        std::cout<<" Done!"<<std::endl;
    }

    if (simData.pzOpts.genVTK) {
        TPZMatModalAnalysis *matPointer =
                dynamic_cast<TPZMatModalAnalysis *>(meshVec[0]->MaterialVec()[1]);
        TPZVec<TPZCompMesh *> temporalMeshVec(2);
        temporalMeshVec[matPointer->H1Index()] = meshVec[1 + matPointer->H1Index()];
        temporalMeshVec[matPointer->HCurlIndex()] =
                meshVec[1 + matPointer->HCurlIndex()];

        std::cout << "Post Processing..." << std::endl;

        TPZStack<std::string> scalnames, vecnames;
        scalnames.Push("Ez");
        vecnames.Push("Et");
        std::string plotfile = simData.pzOpts.prefix + "fieldPlot" + ".vtk";
                                                        // estara na pasta debug
        const int dim = 2;
        an.DefineGraphMesh(dim, scalnames, vecnames,
                           plotfile);  // define malha grafica

        TPZFMatrix<SPZAlwaysComplex<STATE>::type> currentEigenvector(neq,1);
        TPZFMatrix<SPZAlwaysComplex<STATE>::type> scatteredEigen(neqOriginal,1);
        for (int iSol = 0; iSol < eigenValues.size(); iSol++) {
            for(int j = 0; j < eigenVectors.Rows(); j++){
                currentEigenvector(j,0) = simData.pzOpts.absVal ?
                                          std::abs(eigenVectors.GetVal(j,iSol)) :
                                          std::real(eigenVectors.GetVal(j,iSol));
            }
            strmtrx->EquationFilter().Scatter(currentEigenvector, scatteredEigen);
            an.LoadSolution(scatteredEigen);
            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(temporalMeshVec,
                                                               cmesh);
            an.PostProcess(simData.pzOpts.vtkRes);
        }
    }
    gmesh->SetReference(nullptr);
    for (int k = 0; k < meshVec.size(); ++k) {
        meshVec[k]->SetReference(nullptr);
        delete meshVec[k];
        meshVec[k] = nullptr;
    }
    delete gmesh;
    std::cout << "FINISHED!" << std::endl;
    return;
}

void FilterBoundaryEquations(TPZVec<TPZCompMesh *> meshVec,
                             TPZVec<long> &activeEquations, int &neq,
                             int &neqOriginal) {
    TPZCompMesh *cmesh = meshVec[0];

    TPZManVector<long, 1000> allConnects;
    std::set<long> boundConnects;

    for (int iel = 0; iel < cmesh->NElements(); iel++) {
        TPZCompEl *cel = cmesh->ElementVec()[iel];
        if (cel == NULL) {
            continue;
        }
        if (cel->Reference() == NULL) {

            continue;
        }
        TPZBndCond *mat = dynamic_cast<TPZBndCond *>(meshVec[0]->MaterialVec()[cel->Reference()->MaterialId()]);
        if (mat && mat->Type() == 0) {
            std::set<long> boundConnectsEl;
            std::set<long> depBoundConnectsEl;
            std::set<long> indepBoundConnectsEl;
            cel->BuildConnectList(indepBoundConnectsEl, depBoundConnectsEl);
            cel->BuildConnectList(boundConnectsEl);
            for (std::set<long>::iterator iT = boundConnectsEl.begin();
                 iT != boundConnectsEl.end(); iT++) {
                const long val = *iT;
                if (boundConnects.find(val) == boundConnects.end()) {
                    boundConnects.insert(val);
                }
            }
        }
    }
    for (int iCon = 0; iCon < cmesh->NConnects(); iCon++) {
        if (boundConnects.find(iCon) == boundConnects.end()) {
            TPZConnect &con = cmesh->ConnectVec()[iCon];
            if(con.HasDependency())continue;
            int seqnum = con.SequenceNumber();
            int pos = cmesh->Block().Position(seqnum);
            int blocksize = cmesh->Block().Size(seqnum);
            if (blocksize == 0)
                continue;

            int vs = activeEquations.size();
            activeEquations.Resize(vs + blocksize);
            for (int ieq = 0; ieq < blocksize; ieq++) {
                activeEquations[vs + ieq] = pos + ieq;
            }
        }
    }

    neqOriginal = cmesh->NEquations();
    neq = 0;
    TPZCompMesh *cmeshHCurl = meshVec[1];
    TPZCompMesh *cmeshH1 = meshVec[2];
    int nHCurlEquations = 0, nH1Equations = 0;
    TPZMatModalAnalysis *mat =
        dynamic_cast<TPZMatModalAnalysis *>(cmesh->FindMaterial(1));
    for (int iCon = 0; iCon < cmesh->NConnects(); iCon++) {
        bool isH1;
        if (boundConnects.find(iCon) == boundConnects.end()) {
            if (cmesh->ConnectVec()[iCon].HasDependency())
                continue;
            int seqnum = cmesh->ConnectVec()[iCon].SequenceNumber();
            int blocksize = cmesh->Block().Size(seqnum);
            if (mat->H1Index() == 0 && iCon < cmeshH1->NConnects()) {
                isH1 = true;
            } else if (mat->H1Index() == 1 && iCon >= cmeshHCurl->NConnects()) {
                isH1 = true;
            } else {
                isH1 = false;
            }
            for (int ieq = 0; ieq < blocksize; ieq++) {
                neq++;
                isH1 == true ? nH1Equations++ : nHCurlEquations++;
            }
        }
    }
    std::cout << "------\tactive eqs\t-------" << std::endl;
    std::cout << "# H1 equations: " << nH1Equations << std::endl;
    std::cout << "# HCurl equations: " << nHCurlEquations << std::endl;
    std::cout << "# equations: " << neq << std::endl;
    std::cout << "------\t----------\t-------" << std::endl;
    return;
}

//struct used for structural meshing of a general quadrilateral. the coordinates can then be used to create
//triangular or quadrilateral elements.
struct QuadrilateralData {
    const TPZVec<REAL> &coord1;
    const TPZVec<REAL> &coord2;
    const TPZVec<REAL> &coord3;
    const TPZVec<REAL> &coord4;

    std::function<TPZVec<REAL>(const REAL &)> mapSide1;
    std::function<TPZVec<REAL>(const REAL &)> mapSide2;
    std::function<TPZVec<REAL>(const REAL &)> mapSide3;
    std::function<TPZVec<REAL>(const REAL &)> mapSide4;

    bool isSide1nonLinear;
    bool isSide2nonLinear;
    bool isSide3nonLinear;
    bool isSide4nonLinear;

    TPZFMatrix<REAL> side1IntPts;
    TPZFMatrix<REAL> side2IntPts;
    TPZFMatrix<REAL> side3IntPts;
    TPZFMatrix<REAL> side4IntPts;
    TPZFMatrix<REAL> facePts;

    QuadrilateralData(const TPZVec<REAL> &c1,const TPZVec<REAL> &c2,const TPZVec<REAL> &c3,const TPZVec<REAL> &c4) :
            coord1(c1),coord2(c2),coord3(c3),coord4(c4){
        auto map_quad_side_linear = [](const TPZVec<REAL> &coord1 ,const TPZVec<REAL> &coord2,const REAL & s) {
            TPZVec<REAL> point(2,0.);
            point[0] = (coord1[0] - s*coord1[0] + coord2[0] + s*coord2[0])/2.;
            point[1] = (coord1[1] - s*coord1[1] + coord2[1] + s*coord2[1])/2.;
            return point;
        };
        isSide1nonLinear = false;
        mapSide1 = std::bind(map_quad_side_linear,coord1,coord2,std::placeholders::_1);
        isSide2nonLinear = false;
        mapSide2 = std::bind(map_quad_side_linear,coord2,coord3,std::placeholders::_1);
        isSide3nonLinear = false;
        mapSide3 = std::bind(map_quad_side_linear,coord3,coord4,std::placeholders::_1);
        isSide4nonLinear = false;
        mapSide4 = std::bind(map_quad_side_linear,coord4,coord1,std::placeholders::_1);

    }
    void SetMapsLinearity(const bool &m1, const bool &m2, const bool &m3, const bool &m4){
        isSide1nonLinear = m1;
        isSide2nonLinear = m2;
        isSide3nonLinear = m3;
        isSide4nonLinear = m4;
    }

    void CreateQuadrilateral(const int &nQsi, const int &nEta){

        auto getPoints = [this](TPZFMatrix<REAL> &sidePts,const int nPtsQsi,const int nPtsEta,
                                const REAL iniQsi, const REAL endQsi, const REAL iniEta, const REAL endEta){

            sidePts.Resize(nPtsQsi * nPtsEta,2);
            // std::cout<<"nPtsQsi: "<<nPtsQsi <<" nPtsEta : "<< nPtsEta<<std::endl;
            // std::cout<<"iniQsi : "<< iniQsi  <<" endQsi : "<< endQsi<<std::endl;
            // std::cout<<"iniEta : "<< iniEta  <<" endEta : "<< endEta<<std::endl;

            auto mapFace = [this]
                    (const REAL &qsi, const REAL &eta, REAL &x, REAL &y) {
                x=0;
                y=0;
                x += ((1.-qsi)/2)*((1.-eta)/2)*(coord1[0]);
                y += ((1.-qsi)/2)*((1.-eta)/2)*(coord1[1]);
                x += ((1.+qsi)/2)*((1.-eta)/2)*(coord2[0]);
                y += ((1.+qsi)/2)*((1.-eta)/2)*(coord2[1]);
                x += ((1.+qsi)/2)*((1.+eta)/2)*(coord3[0]);
                y += ((1.+qsi)/2)*((1.+eta)/2)*(coord3[1]);
                x += ((1.-qsi)/2)*((1.+eta)/2)*(coord4[0]);
                y += ((1.-qsi)/2)*((1.+eta)/2)*(coord4[1]);
            };

            auto s1 = [](const REAL &qsi, const REAL &eta) { return  1. * qsi;};
            auto s2 = [](const REAL &qsi, const REAL &eta) { return  1. * eta;};
            auto s3 = [](const REAL &qsi, const REAL &eta) { return -1. * qsi;};
            auto s4 = [](const REAL &qsi, const REAL &eta) { return -1. * eta;};

            auto mapLinear = [](const TPZVec<REAL>& coord1side, const TPZVec<REAL>& coord2side,const REAL &s, TPZVec<REAL> &x) {
                x[0] = coord1side[0] * (1.-s)/2 + coord2side[0] * (1+s)/2;
                x[1] = coord1side[1] * (1.-s)/2 + coord2side[1] * (1+s)/2;
            };

            auto correctionFactor1 = [](const REAL &qsi, const REAL &eta) { return 0.25 * (-1 + eta) * (-1 + eta);};
            auto correctionFactor2 = [](const REAL &qsi, const REAL &eta) { return 0.25 * ( 1 + qsi) * ( 1 + qsi);};
            auto correctionFactor3 = [](const REAL &qsi, const REAL &eta) { return 0.25 * ( 1 + eta) * ( 1 + eta);};
            auto correctionFactor4 = [](const REAL &qsi, const REAL &eta) { return 0.25 * (-1 + qsi) * (-1 + qsi);};
            int iPt = 0;
            for(int i = 0; i < nPtsEta ; i++){
                const REAL etaPt = nPtsEta > 1 ? ((REAL)(i) / (nPtsEta - 1)) : 1;
                for(int j = 0; j < nPtsQsi ; j++){
                    const REAL qsiPt = nPtsQsi > 1 ? ((REAL)(j) / (nPtsQsi - 1)) : 1;
                    const REAL qsi = iniQsi * (1. - qsiPt) + endQsi * qsiPt;
                    const REAL eta = iniEta * (1. - etaPt) + endEta * etaPt;
                    REAL xFace;
                    REAL yFace;
                    mapFace(qsi,eta,xFace,yFace);

                    TPZVec<REAL> xSide(2,0.);

                    sidePts(iPt,0) = xFace;
                    sidePts(iPt,1) = yFace;
                    mapLinear(coord1,coord2,s1(qsi,eta),xSide);
                    sidePts(iPt,0) -= isSide1nonLinear ? correctionFactor1(qsi,eta) * (xSide[0] - mapSide1(s1(qsi,eta))[0]) : 0;
                    sidePts(iPt,1) -= isSide1nonLinear ? correctionFactor1(qsi,eta) * (xSide[1] - mapSide1(s1(qsi,eta))[1]) : 0;
                    mapLinear(coord2,coord3,s2(qsi,eta),xSide);
                    sidePts(iPt,0) -= isSide2nonLinear ? correctionFactor2(qsi,eta) * (xSide[0] - mapSide2(s2(qsi,eta))[0]) : 0;
                    sidePts(iPt,1) -= isSide2nonLinear ? correctionFactor2(qsi,eta) * (xSide[1] - mapSide2(s2(qsi,eta))[1]) : 0;
                    mapLinear(coord3,coord4,s3(qsi,eta),xSide);
                    sidePts(iPt,0) -= isSide3nonLinear ? correctionFactor3(qsi,eta) * (xSide[0] - mapSide3(s3(qsi,eta))[0]) : 0;
                    sidePts(iPt,1) -= isSide3nonLinear ? correctionFactor3(qsi,eta) * (xSide[1] - mapSide3(s3(qsi,eta))[1]) : 0;
                    mapLinear(coord4,coord1,s4(qsi,eta),xSide);
                    sidePts(iPt,0) -= isSide4nonLinear ? correctionFactor4(qsi,eta) * (xSide[0] - mapSide4(s4(qsi,eta))[0]) : 0;
                    sidePts(iPt,1) -= isSide4nonLinear ? correctionFactor4(qsi,eta) * (xSide[1] - mapSide4(s4(qsi,eta))[1]) : 0;

                    //std::cout<<"x:   "<<sidePts(iPt,0)<<", y:   "<<sidePts(iPt,1)<<std::endl;
                    iPt++;
                }
            }
        };


        //std::cout<<"--------------- intP ---------------"<<std::endl;
        getPoints(facePts,nQsi,nEta,-1., 1.,-1., 1.);
        const REAL deltaQsiInt = 1./(nQsi - 1);
        const REAL deltaEtaInt = 1./(nEta - 1);
        //std::cout<<"---------------side 1---------------"<<std::endl;
        isSide1nonLinear ?
        getPoints(side1IntPts,(nQsi-1),1,-1.+deltaQsiInt, 1.-deltaQsiInt,-1.,-1.) : (void)side1IntPts.Resize(0,0);
        //std::cout<<"---------------side 2---------------"<<std::endl;
        isSide2nonLinear ?
        getPoints(side2IntPts,1,(nEta-1), 1., 1.,-1.+deltaEtaInt, 1.-deltaEtaInt): (void)side2IntPts.Resize(0,0);
        //std::cout<<"---------------side 3---------------"<<std::endl;
        isSide3nonLinear ?
        getPoints(side3IntPts,(nQsi-1),1, 1.-deltaQsiInt,-1.+deltaQsiInt, 1., 1.): (void)side3IntPts.Resize(0,0);
        //std::cout<<"---------------side 4---------------"<<std::endl;
        isSide4nonLinear ?
        getPoints(side4IntPts,1,(nEta-1),-1.,-1., 1.-deltaEtaInt,-1.+deltaEtaInt): (void)side4IntPts.Resize(0,0);

    }
};


//CreateStepFiberMesh(TPZGeoMesh *&gmesh, const std::string mshFileName, TPZVec<int> &matIdVec, const std::string &prefix,
//            const bool &print, const REAL &scale, const int &factor)
void CreateGMeshRectangularWaveguide(TPZGeoMesh *&gmesh, const std::string mshFileName, TPZVec<int> &matIdVec,
                                     const std::string &prefix, const bool &print, const REAL &scale, const int &factor,
                                     const REAL wDomain, const REAL hDomain, TPZVec<SPZModalAnalysisData::boundtype> &boundTypeVec,
                                     TPZVec<SPZModalAnalysisData::pmltype> &pmlTypeVec, const bool &usingSymmetry,
                                     const SPZModalAnalysisData::boundtype &symmetryType) {
    TPZManVector<int, 3> nx(3, 0);
    TPZManVector<REAL, 3> llCoord(3, 0.), ulCoord(3, 0.), urCoord(3, 0.),
            lrCoord(3, 0.);
    llCoord[0] = 0.;
    llCoord[1] = 0.;

    ulCoord[0] = 0.;
    ulCoord[1] = hDomain*scale;

    urCoord[0] = usingSymmetry ? (wDomain*scale)/2. : wDomain*scale;
    urCoord[1] = hDomain*scale;

    lrCoord[0] = usingSymmetry ? (wDomain*scale)/2. : wDomain*scale;
    lrCoord[1] = 0.;

    nx[0] = factor;
    nx[1] = factor;
    int numl = 1;
    TPZGenGrid *gengrid = new TPZGenGrid(nx, llCoord, urCoord, numl, 0);
    gengrid->SetElementType(ETriangle);

    gmesh = new TPZGeoMesh();
    const int matId = 1; // define id para um material(formulacao fraca)
    const int bc0 = -1;  // define id para um material(cond contorno PEC)
    const int bc1 = -2;  // define id para um material(cond contorno PMC)
    gengrid->Read(gmesh, matId);

    if(usingSymmetry && symmetryType == SPZModalAnalysisData::PMC){
        gengrid->SetBC(gmesh, ulCoord, llCoord, bc1);
    }else{
        gengrid->SetBC(gmesh, ulCoord, llCoord, bc0);
    }
    gengrid->SetBC(gmesh, urCoord, ulCoord, bc0);
    gengrid->SetBC(gmesh, lrCoord, urCoord, bc0);
    gengrid->SetBC(gmesh, llCoord, lrCoord, bc0);

    gmesh->BuildConnectivity();

    if(usingSymmetry && symmetryType == SPZModalAnalysisData::PMC){
        matIdVec.Resize(3);
    }else{
        matIdVec.Resize(2);
    }
    matIdVec[0] = matId;
    matIdVec[1] = bc0;
    if(usingSymmetry && symmetryType == SPZModalAnalysisData::PMC){
        matIdVec[2] = bc1;
    }
    pmlTypeVec.Resize(0);
    if(usingSymmetry && symmetryType == SPZModalAnalysisData::PMC){
        boundTypeVec.Resize(2);
    }else{
        boundTypeVec.Resize(1);
    }
    boundTypeVec[0]=SPZModalAnalysisData::PEC;
    if(usingSymmetry && symmetryType == SPZModalAnalysisData::PMC){
        boundTypeVec[1]=SPZModalAnalysisData::PMC;
    }
#ifdef PZDEBUG
    TPZCheckGeom * Geometrytest = new TPZCheckGeom(gmesh);
    int isBadMeshQ = Geometrytest->PerformCheck();

    if (isBadMeshQ) {
        DebugStop();
    }
#endif

    if(print){
        std::string meshFileName = prefix + "gmesh";
        const size_t strlen = meshFileName.length();
        meshFileName.append(".vtk");
        std::ofstream outVTK(meshFileName.c_str());
        meshFileName.replace(strlen, 4, ".txt");
        std::ofstream outTXT(meshFileName.c_str());

        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outVTK, true);
        gmesh->Print(outTXT);
        outTXT.close();
        outVTK.close();
    }
    delete gengrid;
}

void
CreateStepFiberMesh(TPZGeoMesh *&gmesh, const std::string mshFileName, TPZVec<int> &matIdVec, const std::string &prefix,
                    const bool &print, const REAL &scale, const int &factor,
                    TPZVec<SPZModalAnalysisData::boundtype> &boundTypeVec,
                    TPZVec<SPZModalAnalysisData::pmltype> &pmlTypeVec, const REAL &realRCore, const REAL &dPML,
                    const REAL &boundDist, const REAL &outerReffIndex) {
    const int nDivTCore = factor * 4, nDivRCore = factor * 7, nDivTCladding = factor * 4, nDivPml = factor * 5 + 1;
    //const int nDivTCore = factor * 2, nDivRCore = factor * 3, nDivTCladding = factor * 2, nDivPml = factor * (1+1);
    const int nQuads = 17;
    const int nEdges = 40;



    ///creates refpattern
    TPZAutoPointer<TPZRefPattern> refp;
    char buf[] =
            "4 3 "
            "-50 Quad0000111111111 "
            "-1 -1 0 "
            "1 -1 0 "
            "1 1 0 "
            "-1 1 0 "
            "3 4 0  1  2  3 "
            "2 3 0  1  2 "
            "2 3 0  2  3 ";
    std::istringstream str(buf);
    refp = new TPZRefPattern(str);
    refp->GenerateSideRefPatterns();
    gRefDBase.InsertRefPattern(refp);
    if(!refp)
    {
        DebugStop();
    }
    if(std::min<int>({nDivTCore,nDivRCore,nDivTCladding,nDivPml}) < 2 ) {
        std::cout<<"Mesh has not sufficient divisions."<<std::endl;
        std::cout<<"Minimum is 3."<<std::endl;
        DebugStop();
    }

    const int matIdCore = 1, matIdCladding = 2, matIdBC= 18;
    const int matIdPMLxp = 10,
            matIdPMLyp = 11,
            matIdPMLxm = 12,
            matIdPMLym = 13,
            matIdPMLxpym = 14,
            matIdPMLxpyp = 15,
            matIdPMLxmyp = 16,
            matIdPMLxmym = 17;

    const REAL rCore = scale * realRCore;//8e-6;
    const REAL lengthPML = dPML * outerReffIndex * 2 * M_PI;
    TPZManVector<REAL,2> xc(2, 0.);
    xc[0] = 0.;
    xc[1] = 0.;
    const REAL bound = rCore + boundDist * outerReffIndex * 2 * M_PI;
    TPZManVector<REAL,2> ptCircle_1(2,0.);
    ptCircle_1[0] = xc[0] + rCore * cos(1*M_PI/4);
    ptCircle_1[1] = xc[1] + rCore * sin(1*M_PI/4);
    TPZManVector<REAL,2> ptCircle_2(2,0.);
    ptCircle_2[0] = xc[0] + rCore * cos(3*M_PI/4);
    ptCircle_2[1] = xc[1] + rCore * sin(3*M_PI/4);
    TPZManVector<REAL,2> ptCircle_3(2,0.);
    ptCircle_3[0] = xc[0] + rCore * cos(5*M_PI/4);
    ptCircle_3[1] = xc[1] + rCore * sin(5*M_PI/4);
    TPZManVector<REAL,2> ptCircle_4(2,0.);
    ptCircle_4[0] = xc[0] + rCore * cos(7*M_PI/4);
    ptCircle_4[1] = xc[1] + rCore * sin(7*M_PI/4);

    TPZManVector<REAL,2> ptSquare_1(2,0.);
    ptSquare_1[0] = ptCircle_1[0] / 2.;
    ptSquare_1[1] = ptCircle_1[1] / 2.;
    TPZManVector<REAL,2> ptSquare_2(2,0.);
    ptSquare_2[0] = ptCircle_2[0] / 2.;
    ptSquare_2[1] = ptCircle_2[1] / 2.;
    TPZManVector<REAL,2> ptSquare_3(2,0.);
    ptSquare_3[0] = ptCircle_3[0] / 2.;
    ptSquare_3[1] = ptCircle_3[1] / 2.;
    TPZManVector<REAL,2> ptSquare_4(2,0.);
    ptSquare_4[0] = ptCircle_4[0] / 2.;
    ptSquare_4[1] = ptCircle_4[1] / 2.;

    TPZManVector<REAL,2> ptBound_1(2,0.);
    ptBound_1[0] = 1 * bound;
    ptBound_1[1] = 1 * bound;
    TPZManVector<REAL,2> ptBound_2(2,0.);
    ptBound_2[0] = -1 * bound;
    ptBound_2[1] = 1 * bound;
    TPZManVector<REAL,2> ptBound_3(2,0.);
    ptBound_3[0] = -1 * bound;
    ptBound_3[1] = -1 * bound;
    TPZManVector<REAL,2> ptBound_4(2,0.);
    ptBound_4[0] = 1 * bound;
    ptBound_4[1] = -1 * bound;

    TPZManVector<REAL,2> ptPMLxp_u(2,0.);
    ptPMLxp_u[0] =  1 * bound + lengthPML;
    ptPMLxp_u[1] =  1 * bound;
    TPZManVector<REAL,2> ptPMLxp_l(2,0.);
    ptPMLxp_l[0] =  1 * bound + lengthPML;
    ptPMLxp_l[1] = -1 * bound;

    TPZManVector<REAL,2> ptPMLxm_u(2,0.);
    ptPMLxm_u[0] = -1.*ptPMLxp_u[0];
    ptPMLxm_u[1] = ptPMLxp_u[1];
    TPZManVector<REAL,2> ptPMLxm_l(2,0.);
    ptPMLxm_l[0] = -1.*ptPMLxp_l[0];
    ptPMLxm_l[1] = ptPMLxp_l[1];

    TPZManVector<REAL,2> ptPMLyp_r(2,0.);
    ptPMLyp_r[0] = 1 * bound;
    ptPMLyp_r[1] = 1 * bound + lengthPML;
    TPZManVector<REAL,2> ptPMLyp_l(2,0.);
    ptPMLyp_l[0] = -1 * bound;
    ptPMLyp_l[1] =  1 * bound + lengthPML;

    TPZManVector<REAL,2> ptPMLym_r(2,0.);
    ptPMLym_r[0] =  1. * ptPMLyp_r[0];
    ptPMLym_r[1] = -1. * ptPMLyp_r[1];
    TPZManVector<REAL,2> ptPMLym_l(2,0.);
    ptPMLym_l[0] =  1. * ptPMLyp_l[0];
    ptPMLym_l[1] = -1. * ptPMLyp_l[1];


    TPZManVector<REAL,2> ptPMLxpyp(2,0.);
    ptPMLxpyp[0] = (bound + lengthPML);
    ptPMLxpyp[1] = (bound + lengthPML);
    TPZManVector<REAL,2> ptPMLxmyp(2,0.);
    ptPMLxmyp[0] =-(bound + lengthPML);
    ptPMLxmyp[1] = (bound + lengthPML);
    TPZManVector<REAL,2> ptPMLxmym(2,0.);
    ptPMLxmym[0] =-(bound + lengthPML);
    ptPMLxmym[1] =-(bound + lengthPML);
    TPZManVector<REAL,2> ptPMLxpym(2,0.);
    ptPMLxpym[0] = (bound + lengthPML);
    ptPMLxpym[1] =-(bound + lengthPML);

    auto map_quad_side_arc = [](const TPZVec<REAL> &theta ,const TPZVec<REAL> &xc, const REAL &r, const REAL & s) {
        TPZVec<REAL> point(2,0.);
        point[0] = xc[0] + r*cos((theta[1] + s*theta[1] + theta[0] - s*theta[0])/2.);
        point[1] = xc[1] + r*sin((theta[1] + s*theta[1] + theta[0] - s*theta[0])/2.);
        return point;
    };


    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    // quadrilateral zero - center square                                                                   //
    // quadrilateral one to four - fiber's core, from the right, counter-clockwise                          //
    // quadrilateral five to eight - fiber's cladding, from the right, counter-clockwise                    //
    // quadrilateral nine to twelve - PML regions (excluding corners), from the right, counter-clockwise    //
    // quadrilateral thirteen to sixteen - PML corners, from the top-right, counter-clockwise               //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////CREATE QUADS///////////////////////////////////////
    TPZVec<QuadrilateralData *> quadVec(nQuads,nullptr);
    TPZManVector<int,nQuads> matIdsQuads(nQuads,-1);
    matIdsQuads[0] = matIdCore;

    matIdsQuads[1] = matIdCore;    matIdsQuads[2] = matIdCore;
    matIdsQuads[3] = matIdCore;    matIdsQuads[4] = matIdCore;

    matIdsQuads[5] = matIdCladding;    matIdsQuads[6] = matIdCladding;
    matIdsQuads[7] = matIdCladding;    matIdsQuads[8] = matIdCladding;

    matIdsQuads[9] = matIdPMLxp;    matIdsQuads[10] = matIdPMLyp;
    matIdsQuads[11] = matIdPMLxm;    matIdsQuads[12] = matIdPMLym;
    matIdsQuads[13] = matIdPMLxpyp;    matIdsQuads[14] = matIdPMLxmyp;
    matIdsQuads[15] = matIdPMLxmym;    matIdsQuads[16] = matIdPMLxpym;

    TPZManVector<int,nQuads> nDivQsi(nQuads,-1);
    TPZManVector<int,nQuads> nDivEta(nQuads,-1);
    nDivQsi[0]=nDivRCore;
    nDivEta[0]=nDivRCore;

    nDivQsi[1]=nDivTCore; nDivEta[1]=nDivRCore;
    nDivQsi[2]=nDivRCore; nDivEta[2]=nDivTCore;
    nDivQsi[3]=nDivTCore; nDivEta[3]=nDivRCore;
    nDivQsi[4]=nDivRCore; nDivEta[4]=nDivTCore;

    nDivQsi[5]=nDivTCladding; nDivEta[5]=nDivRCore;
    nDivQsi[6]=nDivRCore;     nDivEta[6]=nDivTCladding;
    nDivQsi[7]=nDivTCladding; nDivEta[7]=nDivRCore;
    nDivQsi[8]=nDivRCore;     nDivEta[8]=nDivTCladding;

    nDivQsi[9]=nDivPml; nDivEta[9]=nDivRCore;
    nDivQsi[10]=nDivRCore;    nDivEta[10]=nDivPml;
    nDivQsi[11]=nDivPml; nDivEta[11]=nDivRCore;
    nDivQsi[12]=nDivRCore;     nDivEta[12]=nDivPml;

    nDivQsi[13]=nDivPml; nDivEta[13]=nDivPml;
    nDivQsi[14]=nDivPml;     nDivEta[14]=nDivPml;
    nDivQsi[15]=nDivPml; nDivEta[15]=nDivPml;
    nDivQsi[16]=nDivPml;     nDivEta[16]=nDivPml;

    TPZManVector<bool,nQuads> side1NonLinearVec(nQuads,false);
    TPZManVector<bool,nQuads> side2NonLinearVec(nQuads,false);
    TPZManVector<bool,nQuads> side3NonLinearVec(nQuads,false);
    TPZManVector<bool,nQuads> side4NonLinearVec(nQuads,false);
    //quad 0 is all false
    side4NonLinearVec[1] = true;
    side1NonLinearVec[2] = true;
    side2NonLinearVec[3] = true;
    side3NonLinearVec[4] = true;

    side2NonLinearVec[5] = true;
    side3NonLinearVec[6] = true;
    side4NonLinearVec[7] = true;
    side1NonLinearVec[8] = true;
    //PMLS are all false

    quadVec[0] = new QuadrilateralData(ptSquare_1,ptSquare_2,ptSquare_3,ptSquare_4);
    quadVec[0]->CreateQuadrilateral(nDivQsi[0],nDivEta[0]);

    TPZVec<REAL> theta(2,0.);

    theta[0] =-1 * M_PI / 4;
    theta[1] = 1 * M_PI / 4;
    quadVec[1] = new QuadrilateralData(ptCircle_1,ptSquare_1,ptSquare_4,ptCircle_4);
    quadVec[1]->mapSide4 = std::bind(map_quad_side_arc,theta, xc, rCore,std::placeholders::_1);
    quadVec[1]->SetMapsLinearity(side1NonLinearVec[1],side2NonLinearVec[1],side3NonLinearVec[1],side4NonLinearVec[1]);
    quadVec[1]->CreateQuadrilateral(nDivQsi[1],nDivEta[1]);

    theta[0] = 1 * M_PI / 4;
    theta[1] = 3 * M_PI / 4;
    quadVec[2] = new QuadrilateralData(ptCircle_1,ptCircle_2,ptSquare_2,ptSquare_1);
    quadVec[2]->mapSide1 = std::bind(map_quad_side_arc,theta, xc, rCore,std::placeholders::_1);
    quadVec[2]->SetMapsLinearity(side1NonLinearVec[2],side2NonLinearVec[2],side3NonLinearVec[2],side4NonLinearVec[2]);
    quadVec[2]->CreateQuadrilateral(nDivQsi[2],nDivEta[2]);

    theta[0] = 3 * M_PI / 4;
    theta[1] = 5 * M_PI / 4;
    quadVec[3] = new QuadrilateralData(ptSquare_2,ptCircle_2,ptCircle_3,ptSquare_3);
    quadVec[3]->mapSide2 = std::bind(map_quad_side_arc,theta, xc, rCore,std::placeholders::_1);
    quadVec[3]->SetMapsLinearity(side1NonLinearVec[3],side2NonLinearVec[3],side3NonLinearVec[3],side4NonLinearVec[3]);
    quadVec[3]->CreateQuadrilateral(nDivQsi[3],nDivEta[3]);

    theta[0] = 5 * M_PI / 4;
    theta[1] = 7 * M_PI / 4;
    quadVec[4] = new QuadrilateralData(ptSquare_4, ptSquare_3,ptCircle_3,ptCircle_4);
    quadVec[4]->mapSide3 = std::bind(map_quad_side_arc,theta, xc, rCore,std::placeholders::_1);
    quadVec[4]->SetMapsLinearity(side1NonLinearVec[4],side2NonLinearVec[4],side3NonLinearVec[4],side4NonLinearVec[4]);
    quadVec[4]->CreateQuadrilateral(nDivQsi[4],nDivEta[4]);

    theta[0] = 1 * M_PI / 4;
    theta[1] =-1 * M_PI / 4;
    quadVec[5] = new QuadrilateralData(ptBound_1,ptCircle_1,ptCircle_4,ptBound_4);
    quadVec[5]->mapSide2 = std::bind(map_quad_side_arc,theta, xc, rCore,std::placeholders::_1);
    quadVec[5]->SetMapsLinearity(side1NonLinearVec[5],side2NonLinearVec[5],side3NonLinearVec[5],side4NonLinearVec[5]);
    quadVec[5]->CreateQuadrilateral(nDivQsi[5],nDivEta[5]);

    theta[0] = 3 * M_PI / 4;
    theta[1] = 1 * M_PI / 4;
    quadVec[6] = new QuadrilateralData(ptBound_1,ptBound_2,ptCircle_2,ptCircle_1);
    quadVec[6]->mapSide3 = std::bind(map_quad_side_arc,theta, xc, rCore,std::placeholders::_1);
    quadVec[6]->SetMapsLinearity(side1NonLinearVec[6],side2NonLinearVec[6],side3NonLinearVec[6],side4NonLinearVec[6]);
    quadVec[6]->CreateQuadrilateral(nDivQsi[6],nDivEta[6]);

    theta[0] = 5 * M_PI / 4;
    theta[1] = 3 * M_PI / 4;
    quadVec[7] = new QuadrilateralData(ptCircle_2,ptBound_2,ptBound_3,ptCircle_3);
    quadVec[7]->mapSide4 = std::bind(map_quad_side_arc,theta, xc, rCore,std::placeholders::_1);
    quadVec[7]->SetMapsLinearity(side1NonLinearVec[7],side2NonLinearVec[7],side3NonLinearVec[7],side4NonLinearVec[7]);
    quadVec[7]->CreateQuadrilateral(nDivQsi[7],nDivEta[7]);

    theta[0] = 7 * M_PI / 4;
    theta[1] = 5 * M_PI / 4;
    quadVec[8] = new QuadrilateralData(ptCircle_4, ptCircle_3,ptBound_3,ptBound_4);
    quadVec[8]->mapSide1 = std::bind(map_quad_side_arc,theta, xc, rCore,std::placeholders::_1);
    quadVec[8]->SetMapsLinearity(side1NonLinearVec[8],side2NonLinearVec[8],side3NonLinearVec[8],side4NonLinearVec[8]);
    quadVec[8]->CreateQuadrilateral(nDivQsi[8],nDivEta[8]);


    quadVec[9] = new QuadrilateralData(ptPMLxp_u,ptBound_1,ptBound_4,ptPMLxp_l);
    quadVec[9]->CreateQuadrilateral(nDivQsi[9],nDivEta[9]);

    quadVec[10] = new QuadrilateralData(ptPMLyp_r,ptPMLyp_l,ptBound_2,ptBound_1);
    quadVec[10]->CreateQuadrilateral(nDivQsi[10],nDivEta[10]);

    quadVec[11] = new QuadrilateralData(ptBound_2,ptPMLxm_u,ptPMLxm_l,ptBound_3);
    quadVec[11]->CreateQuadrilateral(nDivQsi[11],nDivEta[11]);

    quadVec[12] = new QuadrilateralData(ptBound_4,ptBound_3,ptPMLym_l,ptPMLym_r);
    quadVec[12]->CreateQuadrilateral(nDivQsi[12],nDivEta[12]);

    quadVec[13] = new QuadrilateralData(ptPMLxpyp,ptPMLyp_r,ptBound_1,ptPMLxp_u);
    quadVec[13]->CreateQuadrilateral(nDivQsi[13],nDivEta[13]);

    quadVec[14] = new QuadrilateralData(ptPMLyp_l,ptPMLxmyp,ptPMLxm_u,ptBound_2);
    quadVec[14]->CreateQuadrilateral(nDivQsi[14],nDivEta[14]);

    quadVec[15] = new QuadrilateralData(ptBound_3,ptPMLxm_l,ptPMLxmym,ptPMLym_l);
    quadVec[15]->CreateQuadrilateral(nDivQsi[15],nDivEta[15]);

    quadVec[16] = new QuadrilateralData(ptPMLxp_l,ptBound_4,ptPMLym_r,ptPMLxpym);
    quadVec[16]->CreateQuadrilateral(nDivQsi[16],nDivEta[16]);
    ////////////////////////////////////////CREATE NODES///////////////////////////////////////
    ////////////////////////////////////////FOR INTERIOR///////////////////////////////////////
    ///////////////////////////////////////////POINTS//////////////////////////////////////////
    long nodeId = 0;
    gmesh = new TPZGeoMesh;
    TPZManVector<TPZManVector<long,20>,nQuads> elNodeFaceVecs(nQuads,TPZManVector<long,20>(0,-1));

    for(int el = 0; el < nQuads ; el ++) {
        const long nNodesOld = gmesh->NodeVec().NElements();
        const long nNodesEl = (nDivQsi[el]-2)*(nDivEta[el]-2);
        gmesh->NodeVec().Resize(nNodesOld + nNodesEl);
        elNodeFaceVecs[el].Resize(nDivEta[el] * nDivQsi[el],-1);
        //interior nodes only
        for (int i = 1; i < nDivEta[el] - 1; i++) {
            for(int j = 1 ; j < nDivQsi[el] - 1; j++){
                TPZManVector<REAL, 3> node(3, 0.);
                node[0] = quadVec[el]->facePts(i * nDivQsi[el] + j,0);
                node[1] = quadVec[el]->facePts(i * nDivQsi[el] + j,1);
                gmesh->NodeVec()[nodeId].Initialize(node, *gmesh);
                elNodeFaceVecs[el][i * nDivQsi[el] + j] = nodeId;
                nodeId++;
            }
        }
    }

    ////////////////////////////////////////CREATE NODES///////////////////////////////////////
    //////////////////////////////////////////FOR EDGE/////////////////////////////////////////
    ///////////////////////////////////////////POINTS//////////////////////////////////////////
    auto sidePos = [](const int side, const int &i, const int &nQsi, const int &nEta, const bool &reverse){
        const int iP = reverse ? (side%2 ? nQsi -1 - i : nEta -1 - i) : i;
        switch(side){
            case 1 : return iP;
            case 2 : return (iP*nQsi + (nQsi-1));
            case 3 : return (nQsi - 1 - iP + nQsi * (nEta - 1));
            case 4 : return (nEta - 1 - iP)*nQsi;
            default: DebugStop(); return -1;
        }
    };

    auto createEdgeNodes = [sidePos,nDivEta,nDivQsi,quadVec,gmesh]
            (TPZManVector<TPZManVector<long,20>,nQuads> &elNodeFaceVecs, const long &el1, const long &el2,const int &side1, const int &side2, long &nodeId,
             TPZVec<long> &side1pts, TPZVec<long> &side2pts, const bool &revEl2){
        const int nPts = side1 % 2 ? nDivQsi[el1] : nDivEta[el1];
        long nNodesOld = gmesh->NodeVec().NElements();
        long nNodesEl = nPts;
        for (int i = 0; i < nPts; i++) {
            const int posEl1 = sidePos(side1,i,nDivQsi[el1],nDivEta[el1],false);
            const int posEl2 = sidePos(side2,i,nDivQsi[el2],nDivEta[el2],revEl2);
            if(elNodeFaceVecs[el1][posEl1] != -1){
                elNodeFaceVecs[el2][posEl2] = elNodeFaceVecs[el1][posEl1];
                continue;
            }
            if(elNodeFaceVecs[el2][posEl2] != -1){
                elNodeFaceVecs[el1][posEl1] = elNodeFaceVecs[el2][posEl2];
                continue;
            }
            TPZManVector<REAL, 3> node(3, 0.);
            node[0] = quadVec[el1]->facePts(posEl1,0);
            node[1] = quadVec[el1]->facePts(posEl1,1);
            long nNodesOld = gmesh->NodeVec().NElements();
            gmesh->NodeVec().Resize(nNodesOld + 1);
            gmesh->NodeVec()[nodeId].Initialize(node, *gmesh);
            elNodeFaceVecs[el1][posEl1] = nodeId;
            elNodeFaceVecs[el2][posEl2] = nodeId;
            nodeId++;
        }
        TPZFMatrix<REAL> *nonLinearPts = nullptr;
        if(quadVec[el1]->isSide1nonLinear && side1 == 1) { nonLinearPts = &quadVec[el1]->side1IntPts;}
        else if(quadVec[el1]->isSide2nonLinear && side1 == 2) { nonLinearPts = &quadVec[el1]->side2IntPts;}
        else if(quadVec[el1]->isSide3nonLinear && side1 == 3) { nonLinearPts = &quadVec[el1]->side3IntPts;}
        else if(quadVec[el1]->isSide4nonLinear && side1 == 4) { nonLinearPts = &quadVec[el1]->side4IntPts;}
        long nNonLinPts = 0;
        if(nonLinearPts){
            nNonLinPts=nonLinearPts->Rows();
            side1pts.Resize(nNonLinPts);
            side2pts.Resize(nNonLinPts);
            nNodesOld = gmesh->NodeVec().NElements();
            nNodesEl = nNonLinPts;
            gmesh->NodeVec().Resize(nNodesOld + nNodesEl);
        }
        for (int i = 0; i < nNonLinPts; i++){
            const int posEl1 = i;
            const int posEl2 = revEl2? nNonLinPts - 1 - i : i;
            TPZManVector<REAL, 3> node(3, 0.);
            node[0] += (*nonLinearPts)(posEl1,0);
            node[1] += (*nonLinearPts)(posEl1,1);
            gmesh->NodeVec()[nodeId].Initialize(node, *gmesh);
            side1pts[posEl1] = nodeId;
            side2pts[posEl2] = nodeId;
            nodeId++;
        }
    };

    TPZManVector<TPZManVector<long,20>,nQuads> elNodeSide1Vecs(nQuads,TPZManVector<long,20>(0,-1));
    TPZManVector<TPZManVector<long,20>,nQuads> elNodeSide2Vecs(nQuads,TPZManVector<long,20>(0,-1));
    TPZManVector<TPZManVector<long,20>,nQuads> elNodeSide3Vecs(nQuads,TPZManVector<long,20>(0,-1));
    TPZManVector<TPZManVector<long,20>,nQuads> elNodeSide4Vecs(nQuads,TPZManVector<long,20>(0,-1));
    TPZManVector<long,nEdges> quad1Vec(nEdges,-1);//need to keep track of edges and quads
    TPZManVector<long,nEdges> side1Vec(nEdges,-1);
    TPZManVector<long,nEdges> revEl2(nEdges,-1);
    {
        TPZManVector<long,nEdges> quad2Vec(nEdges,-1);
        TPZManVector<long,nEdges> side2Vec(nEdges,-1);
        quad1Vec[0] = 0; quad2Vec[0] = 2; side1Vec[0] = 1; side2Vec[0] = 3; revEl2[0] = true;
        quad1Vec[1] = 0; quad2Vec[1] = 3; side1Vec[1] = 2; side2Vec[1] = 4; revEl2[1] = true;
        quad1Vec[2] = 0; quad2Vec[2] = 4; side1Vec[2] = 3; side2Vec[2] = 1; revEl2[2] = true;
        quad1Vec[3] = 0; quad2Vec[3] = 1; side1Vec[3] = 4; side2Vec[3] = 2; revEl2[3] = true;

        quad1Vec[4] = 1; quad2Vec[4] = 2; side1Vec[4] = 1; side2Vec[4] = 4; revEl2[4] = true;
        quad1Vec[5] = 1; quad2Vec[5] = 4; side1Vec[5] = 3; side2Vec[5] = 4; revEl2[5] = true;
        quad1Vec[6] = 1; quad2Vec[6] = 5; side1Vec[6] = 4; side2Vec[6] = 2; revEl2[6] = true;
        quad1Vec[7] = 2; quad2Vec[7] = 6; side1Vec[7] = 1; side2Vec[7] = 3; revEl2[7] = true;
        quad1Vec[8] = 2; quad2Vec[8] = 3; side1Vec[8] = 2; side2Vec[8] = 1; revEl2[8] = true;
        quad1Vec[9] = 3; quad2Vec[9] = 7; side1Vec[9] = 2; side2Vec[9] = 4; revEl2[9] = true;
        quad1Vec[10] = 3; quad2Vec[10] = 4; side1Vec[10] = 3; side2Vec[10] = 2; revEl2[10] = true;
        quad1Vec[11] = 4; quad2Vec[11] = 8; side1Vec[11] = 3; side2Vec[11] = 1; revEl2[11] = true;

        quad1Vec[12] = 5; quad2Vec[12] = 6; side1Vec[12] = 1; side2Vec[12] = 4; revEl2[12] = true;
        quad1Vec[13] = 5; quad2Vec[13] = 8; side1Vec[13] = 3; side2Vec[13] = 4; revEl2[13] = true;
        quad1Vec[14] = 5; quad2Vec[14] = 9; side1Vec[14] = 4; side2Vec[14] = 2; revEl2[14] = true;
        quad1Vec[15] = 6; quad2Vec[15] =10; side1Vec[15] = 1; side2Vec[15] = 3; revEl2[15] = true;
        quad1Vec[16] = 6; quad2Vec[16] = 7; side1Vec[16] = 2; side2Vec[16] = 1; revEl2[16] = true;
        quad1Vec[17] = 7; quad2Vec[17] =11; side1Vec[17] = 2; side2Vec[17] = 4; revEl2[17] = true;
        quad1Vec[18] = 7; quad2Vec[18] = 8; side1Vec[18] = 3; side2Vec[18] = 2; revEl2[18] = true;
        quad1Vec[19] = 8; quad2Vec[19] =12; side1Vec[19] = 3; side2Vec[19] = 1; revEl2[19] = true;

        quad1Vec[20] = 9; quad2Vec[20] =13; side1Vec[20] = 1; side2Vec[20] = 3; revEl2[20] = true;
        quad1Vec[21] = 9; quad2Vec[21] =16; side1Vec[21] = 3; side2Vec[21] = 1; revEl2[21] = true;
        quad1Vec[22] = 9; quad2Vec[22] = 9; side1Vec[22] = 4; side2Vec[22] = 4; revEl2[22] = false;
        quad1Vec[23] =10; quad2Vec[23] =10; side1Vec[23] = 1; side2Vec[23] = 1; revEl2[23] = false;
        quad1Vec[24] =10; quad2Vec[24] =14; side1Vec[24] = 2; side2Vec[24] = 4; revEl2[24] = true;
        quad1Vec[25] =10; quad2Vec[25] =13; side1Vec[25] = 4; side2Vec[25] = 2; revEl2[25] = true;
        quad1Vec[26] =11; quad2Vec[26] =14; side1Vec[26] = 1; side2Vec[26] = 3; revEl2[26] = true;
        quad1Vec[27] =11; quad2Vec[27] =11; side1Vec[27] = 2; side2Vec[27] = 2; revEl2[27] = false;
        quad1Vec[28] =11; quad2Vec[28] =15; side1Vec[28] = 3; side2Vec[28] = 1; revEl2[28] = true;
        quad1Vec[29] =12; quad2Vec[29] =15; side1Vec[29] = 2; side2Vec[29] = 4; revEl2[29] = true;
        quad1Vec[30] =12; quad2Vec[30] =12; side1Vec[30] = 3; side2Vec[30] = 3; revEl2[30] = false;
        quad1Vec[31] =12; quad2Vec[31] =16; side1Vec[31] = 4; side2Vec[31] = 2; revEl2[31] = true;

        quad1Vec[32] =13; quad2Vec[32] =13; side1Vec[32] = 1; side2Vec[32] = 1; revEl2[32] = false;
        quad1Vec[33] =13; quad2Vec[33] =13; side1Vec[33] = 4; side2Vec[33] = 4; revEl2[33] = false;
        quad1Vec[34] =14; quad2Vec[34] =14; side1Vec[34] = 1; side2Vec[34] = 1; revEl2[34] = false;
        quad1Vec[35] =14; quad2Vec[35] =14; side1Vec[35] = 2; side2Vec[35] = 2; revEl2[35] = false;
        quad1Vec[36] =15; quad2Vec[36] =15; side1Vec[36] = 2; side2Vec[36] = 2; revEl2[36] = false;
        quad1Vec[37] =15; quad2Vec[37] =15; side1Vec[37] = 3; side2Vec[37] = 3; revEl2[37] = false;
        quad1Vec[38] =16; quad2Vec[38] =16; side1Vec[38] = 3; side2Vec[38] = 3; revEl2[38] = false;
        quad1Vec[39] =16; quad2Vec[39] =16; side1Vec[39] = 4; side2Vec[39] = 4; revEl2[39] = false;
        for(int i = 0; i < nEdges; i++){
            TPZVec<long> *side1pts = nullptr;
            TPZVec<long> *side2pts = nullptr;
            switch(side1Vec[i]){
                case 1: side1pts = & elNodeSide1Vecs[quad1Vec[i]]; break;
                case 2: side1pts = & elNodeSide2Vecs[quad1Vec[i]]; break;
                case 3: side1pts = & elNodeSide3Vecs[quad1Vec[i]]; break;
                case 4: side1pts = & elNodeSide4Vecs[quad1Vec[i]]; break;
            }

            switch(side2Vec[i]){
                case 1: side2pts = & elNodeSide1Vecs[quad2Vec[i]]; break;
                case 2: side2pts = & elNodeSide2Vecs[quad2Vec[i]]; break;
                case 3: side2pts = & elNodeSide3Vecs[quad2Vec[i]]; break;
                case 4: side2pts = & elNodeSide4Vecs[quad2Vec[i]]; break;
            }

            createEdgeNodes(elNodeFaceVecs,quad1Vec[i], quad2Vec[i],side1Vec[i],side2Vec[i],nodeId,
                            *side1pts, *side2pts, revEl2[i]);
        }
    }
//    for(int i = 0; i < elNodeFaceVecs.size(); i++){
//        std::cout<<"element "<<i<<std::endl;
//        for(int j = 0; j < elNodeFaceVecs[i].size(); j++){
//            std::cout<<elNodeFaceVecs[i][j]<<"\t";
//        }
//        std::cout<<std::endl;
//    }
    for(int quad = 0; quad < nQuads ; quad ++) {
        for (int i = 0; i < nDivEta[quad] - 1; i++) {
            for(int j = 0 ; j < nDivQsi[quad] - 1; j++){

                const int node0=i * nDivQsi[quad] + j;//Lower left vertex
                const int node1=i * nDivQsi[quad] + j + 1;//Lower right vertex
                const int node2=(i+1) * nDivQsi[quad] + j + 1;//Upper right vertex
                const int node3=(i+1) * nDivQsi[quad] + j;//Upper left vertex

                TPZManVector<long, 3> nodesIdVec(4, 0.);

                const int matId = matIdsQuads[quad];
                nodesIdVec[0] = elNodeFaceVecs[quad][node0];
                nodesIdVec[1] = elNodeFaceVecs[quad][node1];
                nodesIdVec[2] = elNodeFaceVecs[quad][node2];
                nodesIdVec[3] = elNodeFaceVecs[quad][node3];
                if(nodesIdVec[0] == -1 || nodesIdVec[1] == -1 || nodesIdVec[2] == -1 || nodesIdVec[3] == -1){
                    DebugStop();
                }
                if((i == 0 && quadVec[quad]->isSide1nonLinear) ||
                   (j == nDivQsi[quad] - 2 && quadVec[quad]->isSide2nonLinear) ||
                   (i == nDivEta[quad] - 2 && quadVec[quad]->isSide3nonLinear) ||
                   (j == 0 && quadVec[quad]->isSide4nonLinear)){
                    TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> >  * quadEl =
                            new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> >  (nodesIdVec,matId,*gmesh);
                    quadEl->SetRefPattern(refp);
                }
                else{
                    TPZGeoElRefPattern< pzgeom::TPZGeoQuad >  * quadEl =
                            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad >  (nodesIdVec,matId,*gmesh);
                    quadEl->SetRefPattern(refp);
                }
            }
        }
    }

    for(int edge = 0; edge < nEdges ; edge ++) {
        const int quad = quad1Vec[edge];
        const int side = side1Vec[edge];
        if(revEl2[edge] == false){//boundary edge
            const int nArcs = side % 2 ? nDivQsi[quad] -1 : nDivEta[quad] - 1;
            for (int i = 0; i < nArcs; i++) {
                TPZManVector<long, 3> nodesIdVec(2, 0.);
                const int vertex1 = sidePos(side,i,nDivQsi[quad],nDivEta[quad],false);
                const int vertex2 = sidePos(side,i + 1,nDivQsi[quad],nDivEta[quad],false);


                nodesIdVec[0] = elNodeFaceVecs[quad][vertex1];
                nodesIdVec[1] = elNodeFaceVecs[quad][vertex2];
                if(nodesIdVec[0] == -1 || nodesIdVec[1] == -1){
                    DebugStop();
                }
                TPZGeoElRefPattern< pzgeom::TPZGeoLinear > *arc =
                        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (nodesIdVec,matIdBC,*gmesh);
            }
            continue;
        }

        TPZVec<long> *intPointsCoord = nullptr;
        switch(side){
            case 1:
                if(quadVec[quad]->isSide1nonLinear == false) continue;
                intPointsCoord = &(elNodeSide1Vecs[quad]);
                break;
            case 2:
                if(quadVec[quad]->isSide2nonLinear == false) continue;
                intPointsCoord = &(elNodeSide2Vecs[quad]);
                break;
            case 3:
                if(quadVec[quad]->isSide3nonLinear == false) continue;
                intPointsCoord = &(elNodeSide3Vecs[quad]);
                break;
            case 4:
                if(quadVec[quad]->isSide4nonLinear == false) continue;
                intPointsCoord = &(elNodeSide4Vecs[quad]);
                break;
        }

        const int nArcs = side % 2 ? nDivQsi[quad] -1 : nDivEta[quad] - 1;
        //auto sidePos = [](const int side, const int &i, const int &nQsi, const int &nEta, const bool &reverse){
        for (int i = 0; i < nArcs; i++) {
            TPZManVector<long, 3> nodesIdVec(3, 0.);
            const int vertex1 = sidePos(side,i,nDivQsi[quad],nDivEta[quad],false);
            const int vertex2 = sidePos(side,i + 1,nDivQsi[quad],nDivEta[quad],false);


            nodesIdVec[0] = elNodeFaceVecs[quad][vertex1];
            nodesIdVec[1] = elNodeFaceVecs[quad][vertex2];
            nodesIdVec[2] =(*intPointsCoord)[i];//mid-point
            if(nodesIdVec[0] == -1 || nodesIdVec[1] == -1 || nodesIdVec[2] == -1){
                DebugStop();
            }
            TPZGeoElRefPattern< pzgeom::TPZArc3D > *arc =
                    new TPZGeoElRefPattern< pzgeom::TPZArc3D > (nodesIdVec,-24,*gmesh);//random material id
        }
    }
    matIdVec.Resize(11);
    matIdVec[0]=matIdCore;
    matIdVec[1]=matIdCladding;
    pmlTypeVec.Resize(8);
    pmlTypeVec[0]=SPZModalAnalysisData::xp;
    matIdVec[2]=matIdPMLxp;
    pmlTypeVec[1]=SPZModalAnalysisData::yp;
    matIdVec[3]=matIdPMLyp;
    pmlTypeVec[2]=SPZModalAnalysisData::xm;
    matIdVec[4]=matIdPMLxm;
    pmlTypeVec[3]=SPZModalAnalysisData::ym;
    matIdVec[5]=matIdPMLym;
    pmlTypeVec[4]=SPZModalAnalysisData::xpym;
    matIdVec[6]=matIdPMLxpym;
    pmlTypeVec[5]=SPZModalAnalysisData::xpyp;
    matIdVec[7]=matIdPMLxpyp;
    pmlTypeVec[6]=SPZModalAnalysisData::xmyp;
    matIdVec[8]=matIdPMLxmyp;
    pmlTypeVec[7]=SPZModalAnalysisData::xmym;
    matIdVec[9]=matIdPMLxmym;

    boundTypeVec.Resize(1);
    boundTypeVec[0] = SPZModalAnalysisData::PEC;
    matIdVec[10]=matIdBC;

    gmesh->BuildConnectivity();

    TPZVec<TPZGeoEl *> sons;
    //refpatquadtri.rpt
    long nel = gmesh->NElements();
    for (long iel = 0; iel < nel; iel++) {
        TPZGeoEl *gelQuad = gmesh->ElementVec()[iel];
        if(gelQuad->Type() == EQuadrilateral){
            gelQuad->Divide(sons); nel = gmesh->NElements(); continue;
        }
    }

    if(print){
        std::string meshFileName = prefix + "gmesh";
        const size_t strlen = meshFileName.length();
        meshFileName.append(".vtk");
        std::ofstream outVTK(meshFileName.c_str());
        meshFileName.replace(strlen, 4, ".txt");
        std::ofstream outTXT(meshFileName.c_str());

        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outVTK, true);
        gmesh->Print(outTXT);
        outTXT.close();
        outVTK.close();
    }

    return;
}

void
CreateHoleyFiberMesh(TPZGeoMesh *&gmesh, const std::string mshFileName, TPZVec<int> &matIdVec, const std::string &prefix,
                    const bool &print, const REAL &scale, const int &factor,
                    TPZVec<SPZModalAnalysisData::boundtype> &boundTypeVec,
                    TPZVec<SPZModalAnalysisData::pmltype> &pmlTypeVec, const REAL &realRCore, const REAL &dPML,
                    const REAL &boundDist, const REAL &outerReffIndex) {
    const int nDivTHole = factor * 4, nDivRHole = factor * 6,
            nDivCladding1 = factor * 4, nDivCladding2 = factor * 4,nDivCladding3 = factor * 4,
                    nDivPml = factor * 5 + 1;
    const int nQuads = 13;
    const int nEdges = 55;
    const int nPoints = 33;

    if(std::min<int>({nDivTHole,nDivRHole,nDivCladding1,nDivPml}) < 2 ) {
        std::cout<<"Mesh has not sufficient divisions."<<std::endl;
        std::cout<<"Minimum is 2."<<std::endl;
        DebugStop();
    }

    const int matIdHole = 1, matIdCladding = 2, matIdBC= 13;
    const int matIdPMLxp = 10,
            matIdPMLyp = 11,
            matIdPMLxpyp = 12;

    const REAL rCore = scale * 2.5e-6;
    const REAL lengthPML = scale * 2.e-6;
    const REAL bound = scale * 15.75e-6;
    const REAL Lambda = scale * 6.75e-6;
    //holes centers
    TPZManVector<REAL,2> xc1(2, 0.);
    TPZManVector<REAL,2> xc2(2, 0.);
    xc1[0] = Lambda * 0.5;
    xc1[1] = Lambda * std::sqrt(3) * 0.5;
    xc2[0] = Lambda;
    xc2[1] = 0;



    ///all points
    TPZManVector<TPZManVector<REAL,2>,33> pointsVec(nPoints,TPZManVector<REAL,2>(2,0.));
    //hole 1
    pointsVec[0][0]= xc1[0] + 0.5 * sqrt(2)*rCore; pointsVec[0][1]= xc1[1] + 0.5 * sqrt(2)*rCore;
    pointsVec[1][0]= xc1[0] - 0.5 * sqrt(2)*rCore; pointsVec[1][1]= xc1[1] + 0.5 * sqrt(2)*rCore;
    pointsVec[2][0]= xc1[0] - 0.5 * sqrt(2)*rCore; pointsVec[2][1]= xc1[1] - 0.5 * sqrt(2)*rCore;
    pointsVec[3][0]= xc1[0] + 0.5 * sqrt(2)*rCore; pointsVec[3][1]= xc1[1] - 0.5 * sqrt(2)*rCore;

    pointsVec[4][0]= xc1[0] + sqrt(2)*rCore; pointsVec[4][1]= xc1[1] + sqrt(2)*rCore;
    pointsVec[5][0]= xc1[0] - sqrt(2)*rCore; pointsVec[5][1]= xc1[1] + sqrt(2)*rCore;
    pointsVec[6][0]= xc1[0] - sqrt(2)*rCore; pointsVec[6][1]= xc1[1] - sqrt(2)*rCore;
    pointsVec[7][0]= xc1[0] + sqrt(2)*rCore; pointsVec[7][1]= xc1[1] - sqrt(2)*rCore;

    //hole 2
    pointsVec[8][0]= xc2[0] + 0.5 * sqrt(2)*rCore; pointsVec[8][1]= xc2[1] + 0.5 * sqrt(2)*rCore;
    pointsVec[9][0]= xc2[0] - 0.5 * sqrt(2)*rCore; pointsVec[9][1]= xc2[1] + 0.5 * sqrt(2)*rCore;
    pointsVec[10][0]= xc2[0] + 0.5 * sqrt(2)*rCore; pointsVec[10][1]= xc2[1];
    pointsVec[11][0]= xc2[0] - 0.5 * sqrt(2)*rCore; pointsVec[11][1]= xc2[1];

    pointsVec[12][0]= xc2[0] + sqrt(2)*rCore; pointsVec[12][1]= xc2[1] + sqrt(2)*rCore;
    pointsVec[13][0]= xc2[0] - sqrt(2)*rCore; pointsVec[13][1]= xc2[1] + sqrt(2)*rCore;
    pointsVec[14][0]= xc2[0] - rCore; pointsVec[14][1]= xc2[1];
    pointsVec[15][0]= xc2[0] - rCore; pointsVec[15][1]= xc2[1];

    //cladding
    pointsVec[16][0]= 0.; pointsVec[16][1]= 0.;
    pointsVec[17][0]= pointsVec[6][0]; pointsVec[17][1]= 0.;
    pointsVec[18][0]= bound; pointsVec[18][1]= 0.;

    pointsVec[19][0]= pointsVec[16][0]; pointsVec[19][1]= pointsVec[12][1];
    pointsVec[20][0]= pointsVec[17][0]; pointsVec[20][1]= pointsVec[12][1];
    pointsVec[21][0]= pointsVec[18][0]; pointsVec[21][1]= pointsVec[12][1];

    pointsVec[22][0]= pointsVec[16][0]; pointsVec[22][1]= pointsVec[6][1];
    pointsVec[23][0]= pointsVec[17][0]; pointsVec[23][1]= pointsVec[6][1];
    pointsVec[24][0]= pointsVec[18][0]; pointsVec[24][1]= pointsVec[6][1];

    pointsVec[25][0]= pointsVec[16][0]; pointsVec[25][1]= pointsVec[4][1];
    pointsVec[26][0]= pointsVec[17][0]; pointsVec[26][1]= pointsVec[4][1];
    pointsVec[27][0]= pointsVec[18][0]; pointsVec[27][1]= pointsVec[4][1];

    pointsVec[28][0]= pointsVec[16][0]; pointsVec[28][1]= bound;
    pointsVec[29][0]= pointsVec[5][0]; pointsVec[29][1]= bound;
    pointsVec[30][0]= pointsVec[4][0]; pointsVec[30][1]= bound;
    pointsVec[31][0]= pointsVec[12][0]; pointsVec[31][1]= bound;
    pointsVec[32][0]= pointsVec[18][0]; pointsVec[32][1]= bound;
    auto map_quad_side_arc = [](const TPZVec<REAL> &theta ,const TPZVec<REAL> &xc, const REAL &r, const REAL & s) {
        TPZVec<REAL> point(2,0.);
        point[0] = xc[0] + r*cos((theta[1] + s*theta[1] + theta[0] - s*theta[0])/2.);
        point[1] = xc[1] + r*sin((theta[1] + s*theta[1] + theta[0] - s*theta[0])/2.);
        return point;
    };

    TPZFMatrix<int> quadPointsVec(nQuads,4);
    quadPointsVec(0,0) = 0; quadPointsVec(0,1) = 1; quadPointsVec(0,2) = 2; quadPointsVec(0,3) = 3;
    quadPointsVec(1,0) = 4; quadPointsVec(1,1) = 5; quadPointsVec(1,2) = 1; quadPointsVec(1,3) = 0;
    quadPointsVec(2,0) = 1; quadPointsVec(2,1) = 5; quadPointsVec(2,2) = 6; quadPointsVec(2,3) = 2;
    quadPointsVec(3,0) = 4; quadPointsVec(3,1) = 0; quadPointsVec(3,2) = 3; quadPointsVec(3,3) = 7;
    quadPointsVec(4,0) = 3; quadPointsVec(4,1) = 2; quadPointsVec(4,2) = 6; quadPointsVec(4,3) = 7;
    quadPointsVec(5,0) = 8; quadPointsVec(5,1) = 9; quadPointsVec(5,2) = 10; quadPointsVec(5,3) = 11;
    quadPointsVec(6,0) = 12; quadPointsVec(6,1) = 13; quadPointsVec(6,2) = 9; quadPointsVec(6,3) = 8;
    quadPointsVec(7,0) = 9; quadPointsVec(7,1) = 13; quadPointsVec(7,2) = 14; quadPointsVec(7,3) = 10;
    quadPointsVec(8,0) = 12; quadPointsVec(8,1) = 8; quadPointsVec(8,2) = 11; quadPointsVec(8,3) = 15;
    quadPointsVec(9,0) = 20; quadPointsVec(9,1) = 19; quadPointsVec(9,2) = 16; quadPointsVec(9,3) = 17;
    quadPointsVec(10,0) = 13; quadPointsVec(10,1) = 20; quadPointsVec(10,2) = 17; quadPointsVec(10,3) = 14;
    quadPointsVec(11,0) = 21; quadPointsVec(11,1) = 12; quadPointsVec(11,2) = 15; quadPointsVec(11,3) = 18;
    quadPointsVec(12,0) = 6; quadPointsVec(12,1) = 22; quadPointsVec(12,2) = 19; quadPointsVec(12,3) = 20;
    quadPointsVec(13,0) = 7; quadPointsVec(13,1) = 6; quadPointsVec(13,2) = 20; quadPointsVec(13,3) = 13;
    quadPointsVec(14,0) = 23; quadPointsVec(14,1) = 7; quadPointsVec(14,2) = 13; quadPointsVec(14,3) = 12;
    quadPointsVec(15,0) = 24; quadPointsVec(15,1) = 23; quadPointsVec(15,2) = 12; quadPointsVec(15,3) = 21;
    quadPointsVec(16,0) = 5; quadPointsVec(16,1) = 25; quadPointsVec(16,2) = 22; quadPointsVec(16,3) = 6;
    quadPointsVec(17,0) = 26; quadPointsVec(17,1) = 4; quadPointsVec(17,2) = 7; quadPointsVec(17,3) = 23;
    quadPointsVec(10,0) = 27; quadPointsVec(10,1) = 26; quadPointsVec(10,2) = 23; quadPointsVec(10,3) = 24;
    quadPointsVec(19,0) = 29; quadPointsVec(19,1) = 28; quadPointsVec(19,2) = 25; quadPointsVec(19,3) = 5;
    quadPointsVec(20,0) = 30; quadPointsVec(20,1) = 29; quadPointsVec(20,2) = 5; quadPointsVec(20,3) = 4;
    quadPointsVec(21,0) = 31; quadPointsVec(21,1) = 30; quadPointsVec(21,2) = 4; quadPointsVec(21,3) = 26;
    quadPointsVec(22,0) = 32; quadPointsVec(22,1) = 31; quadPointsVec(22,2) = 26; quadPointsVec(22,3) = 27;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    // quadrilateral zero to four - first hole                                                              //
    // quadrilateral five to eight - second hole                                                            //
    // quadrilateral nine to twenty-eight - fiber's cladding                                                //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////CREATE QUADS///////////////////////////////////////
    TPZVec<QuadrilateralData *> quadVec(nQuads,nullptr);
    TPZManVector<int,nQuads> matIdsQuads(nQuads,-1);
    matIdsQuads[0] = matIdHole;
    matIdsQuads[1] = matIdHole;    matIdsQuads[2] = matIdHole;
    matIdsQuads[3] = matIdHole;    matIdsQuads[4] = matIdHole;
    matIdsQuads[5] = matIdHole;    matIdsQuads[6] = matIdHole;
    matIdsQuads[7] = matIdHole;    matIdsQuads[8] = matIdHole;

    matIdsQuads[9] = matIdCladding;    matIdsQuads[10] = matIdCladding;
    matIdsQuads[11] = matIdCladding;    matIdsQuads[12] = matIdCladding;
    matIdsQuads[13] = matIdCladding;    matIdsQuads[14] = matIdCladding;
    matIdsQuads[15] = matIdCladding;    matIdsQuads[16] = matIdCladding;
    matIdsQuads[17] = matIdCladding;    matIdsQuads[18] = matIdCladding;
    matIdsQuads[19] = matIdCladding;    matIdsQuads[20] = matIdCladding;
    matIdsQuads[21] = matIdCladding;    matIdsQuads[22] = matIdCladding;

    TPZManVector<int,nQuads> nDivQsi(nQuads,-1);
    TPZManVector<int,nQuads> nDivEta(nQuads,-1);
    //first hole
    nDivQsi[0]=nDivRHole;nDivEta[0]=nDivRHole;
    nDivQsi[1]=nDivRHole; nDivEta[1]=nDivTHole;
    nDivQsi[2]=nDivTHole; nDivEta[2]=nDivRHole;
    nDivQsi[3]=nDivTHole; nDivEta[3]=nDivRHole;
    nDivQsi[4]=nDivRHole; nDivEta[4]=nDivTHole;
    //second hole
    nDivQsi[5]=nDivRHole; nDivEta[5]=nDivRHole/2;
    nDivQsi[6]=nDivRHole; nDivEta[6]=nDivTHole;
    nDivQsi[7]=nDivTHole; nDivEta[7]=nDivRHole/2;
    nDivQsi[8]=nDivRHole; nDivEta[8]=nDivRHole/2;
    //cladding
    nDivQsi[9]= nDivCladding3; nDivEta[9]=nDivRHole/2;
    nDivQsi[10]=nDivRHole; nDivEta[10]=nDivRHole/2;
    nDivQsi[11]=nDivCladding1; nDivEta[11]=nDivRHole/2;

    nDivQsi[12]=nDivCladding3; nDivEta[12]=nDivCladding2;
    nDivQsi[13]=nDivRHole; nDivEta[13]=nDivCladding2;
    nDivQsi[14]=nDivRHole; nDivEta[14]=nDivCladding2;
    nDivQsi[15]=nDivCladding1; nDivEta[15]=nDivCladding2;

    nDivQsi[16]=nDivCladding3; nDivEta[16]=nDivRHole;
    nDivQsi[17]=nDivRHole; nDivEta[17]=nDivRHole;
    nDivQsi[18]=nDivCladding1; nDivEta[18]=nDivRHole;

    nDivQsi[19]=nDivCladding3; nDivEta[19]=nDivCladding1;
    nDivQsi[20]=nDivRHole; nDivEta[20]=nDivCladding1;
    nDivQsi[21]=nDivRHole; nDivEta[21]=nDivCladding1;
    nDivQsi[22]=nDivCladding1; nDivEta[22]=nDivCladding1;
    TPZManVector<bool,nQuads> side1NonLinearVec(nQuads,false);
    TPZManVector<bool,nQuads> side2NonLinearVec(nQuads,false);
    TPZManVector<bool,nQuads> side3NonLinearVec(nQuads,false);
    TPZManVector<bool,nQuads> side4NonLinearVec(nQuads,false);

    side1NonLinearVec[1] = true;
    side2NonLinearVec[2] = true;
    side4NonLinearVec[3] = true;
    side1NonLinearVec[4] = true;

    side1NonLinearVec[6] = true;
    side2NonLinearVec[7] = true;
    side4NonLinearVec[8] = true;

    side4NonLinearVec[10] = true;
    side2NonLinearVec[11] = true;
    side1NonLinearVec[13] = true;
    side3NonLinearVec[14] = true;
    side4NonLinearVec[16] = true;
    side2NonLinearVec[17] = true;
    side3NonLinearVec[20] = true;

    {
        TPZManVector<TPZManVector<REAL,2>,14> thetaVec(14,TPZManVector<REAL,2>(2,0.));

        thetaVec[0][0] = 1 * M_PI/4;thetaVec[0][1] = 3 * M_PI/4;//quad1
        thetaVec[1][0] = 3 * M_PI/4;thetaVec[1][1] = 5 * M_PI/4;//quad2
        thetaVec[2][0] = 7 * M_PI/4;thetaVec[2][1] = 9 * M_PI/4;//quad3
        thetaVec[3][0] = 5 * M_PI/4;thetaVec[3][1] = 7 * M_PI/4;//quad4

        thetaVec[4][0] = 1 * M_PI/4;thetaVec[4][1] = 3 * M_PI/4;//quad6
        thetaVec[5][0] = 3 * M_PI/4;thetaVec[5][1] = 4 * M_PI/4;//quad7
        thetaVec[6][0] = 0 * M_PI/4;thetaVec[6][1] = 1 * M_PI/4;//quad8

        thetaVec[7][0] = 4 * M_PI/4;thetaVec[7][1] = 3 * M_PI/4;//quad10
        thetaVec[8][0] = 1 * M_PI/4;thetaVec[8][1] = 0 * M_PI/4;//quad11
        thetaVec[9][0] = 7 * M_PI/4;thetaVec[9][1] = 5 * M_PI/4;//quad13
        thetaVec[10][0] = 3 * M_PI/4;thetaVec[10][1] = 1 * M_PI/4;//quad14
        thetaVec[11][0] = 5 * M_PI/4;thetaVec[11][1] = 3 * M_PI/4;//quad16
        thetaVec[12][0] = 9 * M_PI/4;thetaVec[12][1] = 7 * M_PI/4;//quad17
        thetaVec[13][0] = 3 * M_PI/4;thetaVec[13][1] = 1 * M_PI/4;//quad20


        TPZVec<TPZVec<REAL>*> xcRef(14,&xc1);
        xcRef[4] = &xc2;
        xcRef[5] = &xc2;
        xcRef[6] = &xc2;
        xcRef[7] = &xc2;
        xcRef[8] = &xc2;
        xcRef[10] = &xc2;
        int iNonLinear = 0;
        for(int iQuad = 0; iQuad < nQuads; iQuad++){
            quadVec[iQuad] = new QuadrilateralData(
                    pointsVec[quadPointsVec(iQuad,0)],pointsVec[quadPointsVec(iQuad,1)],
                    pointsVec[quadPointsVec(iQuad,2)],pointsVec[quadPointsVec(iQuad,3)]);
            if(side1NonLinearVec[iQuad]){
                TPZVec<REAL> &xc = *(xcRef[iNonLinear]);
                quadVec[iQuad]->mapSide1 = std::bind(map_quad_side_arc,thetaVec[iNonLinear], xc, rCore,std::placeholders::_1);
                iNonLinear++;
            } else if(side2NonLinearVec[iQuad]){
                TPZVec<REAL> &xc = *(xcRef[iNonLinear]);
                quadVec[iQuad]->mapSide2 = std::bind(map_quad_side_arc,thetaVec[iNonLinear], xc, rCore,std::placeholders::_1);
                iNonLinear++;
            } else if(side3NonLinearVec[iQuad]){
                TPZVec<REAL> &xc = *(xcRef[iNonLinear]);
                quadVec[iQuad]->mapSide3 = std::bind(map_quad_side_arc,thetaVec[iNonLinear], xc, rCore,std::placeholders::_1);
                iNonLinear++;
            } else if(side4NonLinearVec[iQuad]){
                TPZVec<REAL> &xc = *(xcRef[iNonLinear]);
                quadVec[iQuad]->mapSide4 = std::bind(map_quad_side_arc,thetaVec[iNonLinear], xc, rCore,std::placeholders::_1);
                iNonLinear++;
            }
            quadVec[iQuad]->SetMapsLinearity(side1NonLinearVec[iQuad],side2NonLinearVec[iQuad],side3NonLinearVec[iQuad],side4NonLinearVec[iQuad]);
            quadVec[iQuad]->CreateQuadrilateral(nDivQsi[iQuad],nDivEta[iQuad]);
        }
    }

    ////////////////////////////////////////CREATE NODES///////////////////////////////////////
    ////////////////////////////////////////FOR INTERIOR///////////////////////////////////////
    ///////////////////////////////////////////POINTS//////////////////////////////////////////
    long nodeId = 0;
    gmesh = new TPZGeoMesh;
    TPZManVector<TPZManVector<long,20>,nQuads> elNodeFaceVecs(nQuads,TPZManVector<long,20>(0,-1));

    for(int el = 0; el < nQuads ; el ++) {
        const long nNodesOld = gmesh->NodeVec().NElements();
        const long nNodesEl = (nDivQsi[el]-2)*(nDivEta[el]-2);
        gmesh->NodeVec().Resize(nNodesOld + nNodesEl);
        elNodeFaceVecs[el].Resize(nDivEta[el] * nDivQsi[el],-1);
        //interior nodes only
        for (int i = 1; i < nDivEta[el] - 1; i++) {
            for(int j = 1 ; j < nDivQsi[el] - 1; j++){
                TPZManVector<REAL, 3> node(3, 0.);
                node[0] = quadVec[el]->facePts(i * nDivQsi[el] + j,0);
                node[1] = quadVec[el]->facePts(i * nDivQsi[el] + j,1);
                gmesh->NodeVec()[nodeId].Initialize(node, *gmesh);
                elNodeFaceVecs[el][i * nDivQsi[el] + j] = nodeId;
                nodeId++;
            }
        }
    }

    ////////////////////////////////////////CREATE NODES///////////////////////////////////////
    //////////////////////////////////////////FOR EDGE/////////////////////////////////////////
    ///////////////////////////////////////////POINTS//////////////////////////////////////////
    auto sidePos = [](const int side, const int &i, const int &nQsi, const int &nEta, const bool &reverse){
        const int iP = reverse ? (side%2 ? nQsi -1 - i : nEta -1 - i) : i;
        switch(side){
            case 1 : return iP;
            case 2 : return (iP*nQsi + (nQsi-1));
            case 3 : return (nQsi - 1 - iP + nQsi * (nEta - 1));
            case 4 : return (nEta - 1 - iP)*nQsi;
            default: DebugStop(); return -1;
        }
    };

    auto createEdgeNodes = [sidePos,nDivEta,nDivQsi,quadVec,gmesh]
            (TPZManVector<TPZManVector<long,20>,nQuads> &elNodeFaceVecs, const long &el1, const long &el2,const int &side1, const int &side2, long &nodeId,
             TPZVec<long> &side1pts, TPZVec<long> &side2pts, const bool &revEl2){
        const int nPts = side1 % 2 ? nDivQsi[el1] : nDivEta[el1];
        long nNodesOld = gmesh->NodeVec().NElements();
        long nNodesEl = nPts;
        for (int i = 0; i < nPts; i++) {
            const int posEl1 = sidePos(side1,i,nDivQsi[el1],nDivEta[el1],false);
            const int posEl2 = sidePos(side2,i,nDivQsi[el2],nDivEta[el2],revEl2);
            if(elNodeFaceVecs[el1][posEl1] != -1){
                elNodeFaceVecs[el2][posEl2] = elNodeFaceVecs[el1][posEl1];
                continue;
            }
            if(elNodeFaceVecs[el2][posEl2] != -1){
                elNodeFaceVecs[el1][posEl1] = elNodeFaceVecs[el2][posEl2];
                continue;
            }
            TPZManVector<REAL, 3> node(3, 0.);
            node[0] = quadVec[el1]->facePts(posEl1,0);
            node[1] = quadVec[el1]->facePts(posEl1,1);
            long nNodesOld = gmesh->NodeVec().NElements();
            gmesh->NodeVec().Resize(nNodesOld + 1);
            gmesh->NodeVec()[nodeId].Initialize(node, *gmesh);
            elNodeFaceVecs[el1][posEl1] = nodeId;
            elNodeFaceVecs[el2][posEl2] = nodeId;
            nodeId++;
        }
        TPZFMatrix<REAL> *nonLinearPts = nullptr;
        if(quadVec[el1]->isSide1nonLinear && side1 == 1) { nonLinearPts = &quadVec[el1]->side1IntPts;}
        else if(quadVec[el1]->isSide2nonLinear && side1 == 2) { nonLinearPts = &quadVec[el1]->side2IntPts;}
        else if(quadVec[el1]->isSide3nonLinear && side1 == 3) { nonLinearPts = &quadVec[el1]->side3IntPts;}
        else if(quadVec[el1]->isSide4nonLinear && side1 == 4) { nonLinearPts = &quadVec[el1]->side4IntPts;}
        long nNonLinPts = 0;
        if(nonLinearPts){
            nNonLinPts=nonLinearPts->Rows();
            side1pts.Resize(nNonLinPts);
            side2pts.Resize(nNonLinPts);
            nNodesOld = gmesh->NodeVec().NElements();
            nNodesEl = nNonLinPts;
            gmesh->NodeVec().Resize(nNodesOld + nNodesEl);
        }
        for (int i = 0; i < nNonLinPts; i++){
            const int posEl1 = i;
            const int posEl2 = revEl2? nNonLinPts - 1 - i : i;
            TPZManVector<REAL, 3> node(3, 0.);
            node[0] += (*nonLinearPts)(posEl1,0);
            node[1] += (*nonLinearPts)(posEl1,1);
            gmesh->NodeVec()[nodeId].Initialize(node, *gmesh);
            side1pts[posEl1] = nodeId;
            side2pts[posEl2] = nodeId;
            nodeId++;
        }
    };

    TPZManVector<TPZManVector<long,20>,nQuads> elNodeSide1Vecs(nQuads,TPZManVector<long,20>(0,-1));
    TPZManVector<TPZManVector<long,20>,nQuads> elNodeSide2Vecs(nQuads,TPZManVector<long,20>(0,-1));
    TPZManVector<TPZManVector<long,20>,nQuads> elNodeSide3Vecs(nQuads,TPZManVector<long,20>(0,-1));
    TPZManVector<TPZManVector<long,20>,nQuads> elNodeSide4Vecs(nQuads,TPZManVector<long,20>(0,-1));
    TPZManVector<long,nEdges> quad1Vec(nEdges,-1);//need to keep track of edges and quads
    TPZManVector<long,nEdges> side1Vec(nEdges,-1);
    TPZManVector<long,nEdges> revEl2(nEdges,-1);
    {
        TPZManVector<long,nEdges> quad2Vec(nEdges,-1);
        TPZManVector<long,nEdges> side2Vec(nEdges,-1);
      //quad1Vec[0] = ok; quad2Vec[0] = ok; side1Vec[0] =nok; side2Vec[0] =nok; revEl2[0] = ok;
        quad1Vec[0] = 0; quad2Vec[0] = 1; side1Vec[0] = 1; side2Vec[0] = 3; revEl2[0] = true;
        quad1Vec[1] = 0; quad2Vec[1] = 2; side1Vec[1] = 2; side2Vec[1] = 4; revEl2[1] = true;
        quad1Vec[2] = 0; quad2Vec[2] = 4; side1Vec[2] = 3; side2Vec[2] = 1; revEl2[2] = true;
        quad1Vec[3] = 0; quad2Vec[3] = 3; side1Vec[3] = 4; side2Vec[3] = 2; revEl2[3] = true;

        quad1Vec[4] = 1; quad2Vec[4] = 3; side1Vec[4] = 1; side2Vec[4] = 4; revEl2[4] = true;
        quad1Vec[5] = 1; quad2Vec[5] = 2; side1Vec[5] = 3; side2Vec[5] = 4; revEl2[5] = true;

        quad1Vec[6] = 2; quad2Vec[6] = 4; side1Vec[6] = 4; side2Vec[6] = 2; revEl2[6] = true;

        quad1Vec[7] = 4; quad2Vec[7] = 3; side1Vec[7] = 1; side2Vec[7] = 3; revEl2[7] = true;

        quad1Vec[8] = 1; quad2Vec[8] =20; side1Vec[8] = 2; side2Vec[8] = 1; revEl2[8] = true;
        quad1Vec[9] = 2; quad2Vec[9] =16; side1Vec[9] = 2; side2Vec[9] = 4; revEl2[9] = true;
        quad1Vec[10] = 4; quad2Vec[10] =13; side1Vec[10] = 3; side2Vec[10] = 2; revEl2[10] = true;
        quad1Vec[11] = 3; quad2Vec[11] =17; side1Vec[11] = 3; side2Vec[11] = 1; revEl2[11] = true;

        quad1Vec[12] = 5; quad2Vec[12] = 6; side1Vec[12] = 1; side2Vec[12] = 4; revEl2[12] = true;
        quad1Vec[13] = 5; quad2Vec[13] = 7; side1Vec[13] = 3; side2Vec[13] = 4; revEl2[13] = true;
        quad1Vec[14] = 5; quad2Vec[14] = 8; side1Vec[14] = 4; side2Vec[14] = 2; revEl2[14] = true;

        quad1Vec[15] = 6; quad2Vec[15] = 8; side1Vec[15] = 1; side2Vec[15] = 3; revEl2[15] = true;
        quad1Vec[16] = 6; quad2Vec[16] = 7; side1Vec[16] = 2; side2Vec[16] = 1; revEl2[16] = true;
        quad1Vec[17] = 6; quad2Vec[17] =14; side1Vec[17] = 2; side2Vec[17] = 4; revEl2[17] = true;
        quad1Vec[18] = 7; quad2Vec[18] =10; side1Vec[18] = 3; side2Vec[18] = 2; revEl2[18] = true;
        quad1Vec[19] = 8; quad2Vec[19] =11; side1Vec[19] = 3; side2Vec[19] = 1; revEl2[19] = true;
        quad1Vec[20] = 9; quad2Vec[20] =12; side1Vec[20] = 1; side2Vec[20] = 3; revEl2[20] = true;
        quad1Vec[21] = 9; quad2Vec[21] =10; side1Vec[21] = 3; side2Vec[21] = 1; revEl2[21] = true;
        quad1Vec[22] =10; quad2Vec[22] =13; side1Vec[22] = 4; side2Vec[22] = 4; revEl2[22] = true;
        quad1Vec[23] =11; quad2Vec[23] =15; side1Vec[23] = 1; side2Vec[23] = 1; revEl2[23] = true;
        quad1Vec[24] =12; quad2Vec[24] =16; side1Vec[24] = 2; side2Vec[24] = 4; revEl2[24] = true;
        quad1Vec[25] =12; quad2Vec[25] =13; side1Vec[25] = 4; side2Vec[25] = 2; revEl2[25] = true;
        quad1Vec[26] =13; quad2Vec[26] =14; side1Vec[26] = 1; side2Vec[26] = 3; revEl2[26] = true;
        quad1Vec[27] =14; quad2Vec[27] =17; side1Vec[27] = 2; side2Vec[27] = 2; revEl2[27] = true;
        quad1Vec[28] =14; quad2Vec[28] =15; side1Vec[28] = 3; side2Vec[28] = 1; revEl2[28] = true;
        quad1Vec[29] =15; quad2Vec[29] =18; side1Vec[29] = 2; side2Vec[29] = 4; revEl2[29] = true;
        quad1Vec[30] =16; quad2Vec[30] =19; side1Vec[30] = 3; side2Vec[30] = 3; revEl2[30] = true;
        quad1Vec[31] =17; quad2Vec[31] =21; side1Vec[31] = 4; side2Vec[31] = 2; revEl2[31] = true;
        quad1Vec[32] =17; quad2Vec[32] =18; side1Vec[32] = 1; side2Vec[32] = 1; revEl2[32] = true;
        quad1Vec[33] =18; quad2Vec[33] =22; side1Vec[33] = 4; side2Vec[33] = 4; revEl2[33] = true;
        quad1Vec[34] =19; quad2Vec[34] =20; side1Vec[34] = 1; side2Vec[34] = 1; revEl2[34] = true;
        quad1Vec[35] =20; quad2Vec[35] =21; side1Vec[35] = 2; side2Vec[35] = 2; revEl2[35] = true;
        quad1Vec[36] =21; quad2Vec[36] =22; side1Vec[36] = 2; side2Vec[36] = 2; revEl2[36] = true;


        quad1Vec[37] =19; quad2Vec[37] =19; side1Vec[37] = 3; side2Vec[37] = 3; revEl2[37] = false;
        quad1Vec[38] =16; quad2Vec[38] =16; side1Vec[38] = 3; side2Vec[38] = 3; revEl2[38] = false;
        quad1Vec[39] =12; quad2Vec[39] =12; side1Vec[39] = 4; side2Vec[39] = 4; revEl2[39] = false;
        quad1Vec[40] = 9; quad2Vec[40] = 9; side1Vec[40] = 4; side2Vec[40] = 4; revEl2[40] = false;

        quad1Vec[41] = 9; quad2Vec[41] = 9; side1Vec[41] = 3; side2Vec[41] = 3; revEl2[41] = false;
        quad1Vec[42] =10; quad2Vec[42] =10; side1Vec[42] = 3; side2Vec[42] = 3; revEl2[42] = false;
        quad1Vec[43] = 7; quad2Vec[43] = 7; side1Vec[43] = 4; side2Vec[43] = 4; revEl2[43] = false;
        quad1Vec[44] = 5; quad2Vec[44] = 5; side1Vec[44] = 4; side2Vec[44] = 4; revEl2[44] = false;
        quad1Vec[45] = 8; quad2Vec[45] = 8; side1Vec[45] = 4; side2Vec[45] = 4; revEl2[45] = false;
        quad1Vec[46] =11; quad2Vec[46] =11; side1Vec[46] = 4; side2Vec[46] = 4; revEl2[46] = false;

        quad1Vec[47] =11; quad2Vec[47] =11; side1Vec[47] = 3; side2Vec[47] = 3; revEl2[47] = false;
        quad1Vec[48] =15; quad2Vec[48] =15; side1Vec[48] = 3; side2Vec[48] = 3; revEl2[48] = false;
        quad1Vec[49] =18; quad2Vec[49] =18; side1Vec[49] = 4; side2Vec[49] = 4; revEl2[49] = false;
        quad1Vec[50] =22; quad2Vec[50] =22; side1Vec[50] = 4; side2Vec[50] = 4; revEl2[50] = false;

        quad1Vec[51] =22; quad2Vec[51] =22; side1Vec[51] = 3; side2Vec[51] = 3; revEl2[51] = false;
        quad1Vec[52] =21; quad2Vec[52] =21; side1Vec[52] = 3; side2Vec[52] = 3; revEl2[52] = false;
        quad1Vec[53] =19; quad2Vec[53] =19; side1Vec[53] = 4; side2Vec[53] = 4; revEl2[53] = false;
        quad1Vec[54] =20; quad2Vec[54] =20; side1Vec[54] = 4; side2Vec[54] = 4; revEl2[54] = false;
        for(int i = 0; i < nEdges; i++){
            TPZVec<long> *side1pts = nullptr;
            TPZVec<long> *side2pts = nullptr;
            switch(side1Vec[i]){
                case 1: side1pts = & elNodeSide1Vecs[quad1Vec[i]]; break;
                case 2: side1pts = & elNodeSide2Vecs[quad1Vec[i]]; break;
                case 3: side1pts = & elNodeSide3Vecs[quad1Vec[i]]; break;
                case 4: side1pts = & elNodeSide4Vecs[quad1Vec[i]]; break;
            }

            switch(side2Vec[i]){
                case 1: side2pts = & elNodeSide1Vecs[quad2Vec[i]]; break;
                case 2: side2pts = & elNodeSide2Vecs[quad2Vec[i]]; break;
                case 3: side2pts = & elNodeSide3Vecs[quad2Vec[i]]; break;
                case 4: side2pts = & elNodeSide4Vecs[quad2Vec[i]]; break;
            }

            createEdgeNodes(elNodeFaceVecs,quad1Vec[i], quad2Vec[i],side1Vec[i],side2Vec[i],nodeId,
                            *side1pts, *side2pts, revEl2[i]);
        }
    }
///creates refpattern
    TPZAutoPointer<TPZRefPattern> refp;
    char buf[] =
            "4 3 "
            "-50 Quad0000111111111 "
            "-1 -1 0 "
            "1 -1 0 "
            "1 1 0 "
            "-1 1 0 "
            "3 4 0  1  2  3 "
            "2 3 0  1  2 "
            "2 3 0  2  3 ";
    std::istringstream str(buf);
    refp = new TPZRefPattern(str);
    refp->GenerateSideRefPatterns();
    gRefDBase.InsertRefPattern(refp);
    if(!refp)
    {
        DebugStop();
    }
    for(int quad = 0; quad < nQuads ; quad ++) {
        for (int i = 0; i < nDivEta[quad] - 1; i++) {
            for(int j = 0 ; j < nDivQsi[quad] - 1; j++){

                const int node0=i * nDivQsi[quad] + j;//Lower left vertex
                const int node1=i * nDivQsi[quad] + j + 1;//Lower right vertex
                const int node2=(i+1) * nDivQsi[quad] + j + 1;//Upper right vertex
                const int node3=(i+1) * nDivQsi[quad] + j;//Upper left vertex

                TPZManVector<long, 3> nodesIdVec(4, 0.);

                const int matId = matIdsQuads[quad];
                nodesIdVec[0] = elNodeFaceVecs[quad][node0];
                nodesIdVec[1] = elNodeFaceVecs[quad][node1];
                nodesIdVec[2] = elNodeFaceVecs[quad][node2];
                nodesIdVec[3] = elNodeFaceVecs[quad][node3];
                if(nodesIdVec[0] == -1 || nodesIdVec[1] == -1 || nodesIdVec[2] == -1 || nodesIdVec[3] == -1){
                    DebugStop();
                }
                if((i == 0 && quadVec[quad]->isSide1nonLinear) ||
                   (j == nDivQsi[quad] - 2 && quadVec[quad]->isSide2nonLinear) ||
                   (i == nDivEta[quad] - 2 && quadVec[quad]->isSide3nonLinear) ||
                   (j == 0 && quadVec[quad]->isSide4nonLinear)){
                    TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> >  * quadEl =
                            new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> >  (nodesIdVec,matId,*gmesh);
                    quadEl->SetRefPattern(refp);
                }
                else{
                    TPZGeoElRefPattern< pzgeom::TPZGeoQuad >  * quadEl =
                            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad >  (nodesIdVec,matId,*gmesh);
                    quadEl->SetRefPattern(refp);
                }
            }
        }
    }

    for(int edge = 0; edge < nEdges ; edge ++) {
        const int quad = quad1Vec[edge];
        const int side = side1Vec[edge];
        if(revEl2[edge] == false){//boundary edge
            const int nArcs = side % 2 ? nDivQsi[quad] -1 : nDivEta[quad] - 1;
            for (int i = 0; i < nArcs; i++) {
                TPZManVector<long, 3> nodesIdVec(2, 0.);
                const int vertex1 = sidePos(side,i,nDivQsi[quad],nDivEta[quad],false);
                const int vertex2 = sidePos(side,i + 1,nDivQsi[quad],nDivEta[quad],false);


                nodesIdVec[0] = elNodeFaceVecs[quad][vertex1];
                nodesIdVec[1] = elNodeFaceVecs[quad][vertex2];
                if(nodesIdVec[0] == -1 || nodesIdVec[1] == -1){
                    DebugStop();
                }
                TPZGeoElRefPattern< pzgeom::TPZGeoLinear > *arc =
                        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (nodesIdVec,matIdBC,*gmesh);
            }
            continue;
        }

        TPZVec<long> *intPointsCoord = nullptr;
        switch(side){
            case 1:
                if(quadVec[quad]->isSide1nonLinear == false) continue;
                intPointsCoord = &(elNodeSide1Vecs[quad]);
                break;
            case 2:
                if(quadVec[quad]->isSide2nonLinear == false) continue;
                intPointsCoord = &(elNodeSide2Vecs[quad]);
                break;
            case 3:
                if(quadVec[quad]->isSide3nonLinear == false) continue;
                intPointsCoord = &(elNodeSide3Vecs[quad]);
                break;
            case 4:
                if(quadVec[quad]->isSide4nonLinear == false) continue;
                intPointsCoord = &(elNodeSide4Vecs[quad]);
                break;
        }

        const int nArcs = side % 2 ? nDivQsi[quad] -1 : nDivEta[quad] - 1;
        //auto sidePos = [](const int side, const int &i, const int &nQsi, const int &nEta, const bool &reverse){
        for (int i = 0; i < nArcs; i++) {
            TPZManVector<long, 3> nodesIdVec(3, 0.);
            const int vertex1 = sidePos(side,i,nDivQsi[quad],nDivEta[quad],false);
            const int vertex2 = sidePos(side,i + 1,nDivQsi[quad],nDivEta[quad],false);


            nodesIdVec[0] = elNodeFaceVecs[quad][vertex1];
            nodesIdVec[1] = elNodeFaceVecs[quad][vertex2];
            nodesIdVec[2] =(*intPointsCoord)[i];//mid-point
            if(nodesIdVec[0] == -1 || nodesIdVec[1] == -1 || nodesIdVec[2] == -1){
                DebugStop();
            }
            TPZGeoElRefPattern< pzgeom::TPZArc3D > *arc =
                    new TPZGeoElRefPattern< pzgeom::TPZArc3D > (nodesIdVec,-24,*gmesh);//random material id
        }
    }
    matIdVec.Resize(3);
    matIdVec[0]=matIdHole;
    matIdVec[1]=matIdCladding;
    pmlTypeVec.Resize(0);
//    pmlTypeVec.Resize(3);
//    pmlTypeVec[0]=SPZModalAnalysisData::xp;
//    matIdVec[2]=matIdPMLxp;
//    pmlTypeVec[1]=SPZModalAnalysisData::yp;
//    matIdVec[3]=matIdPMLyp;
//    pmlTypeVec[2]=SPZModalAnalysisData::xpyp;
//    matIdVec[4]=matIdPMLxpyp;

    boundTypeVec.Resize(1);
    boundTypeVec[0] = SPZModalAnalysisData::PEC;
    matIdVec[matIdVec.size()-1]=matIdBC;

    gmesh->BuildConnectivity();

    TPZVec<TPZGeoEl *> sons;
    //refpatquadtri.rpt
    long nel = gmesh->NElements();
    for (long iel = 0; iel < nel; iel++) {
        TPZGeoEl *gelQuad = gmesh->ElementVec()[iel];
        if(gelQuad->Type() == EQuadrilateral){
            gelQuad->Divide(sons); nel = gmesh->NElements(); continue;
        }
    }

    if(print){
        std::string meshFileName = prefix + "gmesh";
        const size_t strlen = meshFileName.length();
        meshFileName.append(".vtk");
        std::ofstream outVTK(meshFileName.c_str());
        meshFileName.replace(strlen, 4, ".txt");
        std::ofstream outTXT(meshFileName.c_str());

        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outVTK, true);
        gmesh->Print(outTXT);
        outTXT.close();
        outVTK.close();
    }

    return;
}
void
ReadGMesh(TPZGeoMesh *&gmesh, const std::string mshFileName, TPZVec<int> &matIdVec, const std::string &prefix,
          const bool &print, const REAL &scale) {
    TPZGmshReader meshReader;
    meshReader.SetfDimensionlessL(scale);
    gmesh = meshReader.GeometricGmshMesh(mshFileName);
#ifdef PZDEBUG
//	TPZCheckGeom * Geometrytest = new TPZCheckGeom(gmesh);
//	int isBadMeshQ = Geometrytest->PerformCheck();
//
//	if (isBadMeshQ) {
//		DebugStop();
//	}
#endif
    auto matIds = meshReader.fMaterialDataVec.fMatID;
    matIdVec.Resize(matIds.size());
    int i = 0;
    for(auto id = matIds.begin(); id != matIds.end(); id++,i++ )    {
      matIdVec[i] = *id;
    }
    if(print){
        std::string meshFileName = prefix + "gmesh";
        const size_t strlen = meshFileName.length();
        meshFileName.append(".vtk");
        std::ofstream outVTK(meshFileName.c_str());
        meshFileName.replace(strlen, 4, ".txt");
        std::ofstream outTXT(meshFileName.c_str());

        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outVTK, true);
        gmesh->Print(outTXT);
        outTXT.close();
        outVTK.close();
    }

    return;
}

void CreateCMesh(TPZVec<TPZCompMesh *> &meshVecOut, TPZGeoMesh *gmesh, int pOrder, TPZVec<int> &matIdVec,
                 TPZVec<STATE> &urVec, TPZVec<STATE> &erVec, REAL f0, bool isCutOff, const std::string &prefix,
                 const bool &print, const REAL &scale, const REAL &alphaMax,
                 TPZVec<SPZModalAnalysisData::boundtype> &boundTypeVec,
                 TPZVec<SPZModalAnalysisData::pmltype> &pmlTypeVec) {
    const int dim = 2;

    TPZManVector<int,8> volMatIdVec(matIdVec.size()-boundTypeVec.size()-pmlTypeVec.size());
    for(int i = 0; i< volMatIdVec.size(); i++){
      volMatIdVec[i] = matIdVec[i];
    }
    TPZManVector<int,8>pmlMatIdVec(pmlTypeVec.size());
    for(int i = 0; i< pmlMatIdVec.size(); i++){
        pmlMatIdVec[i] = matIdVec[volMatIdVec.size()+i];
    }
    TPZManVector<int,8>boundMatIdVec(boundTypeVec.size());
    for(int i = 0; i< boundMatIdVec.size(); i++){
        boundMatIdVec[i] = matIdVec[volMatIdVec.size()+pmlMatIdVec.size()+i];
    }

    if(volMatIdVec.size() !=urVec.size()){
        std::cout<<"The number of materials is not consistent with the mesh materials."<<std::endl;
        std::cout<<"The number of materials is "<<urVec.size()<<" and there are "<<volMatIdVec.size();
        std::cout<<" mesh materials."<<std::endl;
        DebugStop();
    }
    if(pmlMatIdVec.size() != pmlTypeVec.size()){
        std::cout<<"The number of PML types is not consistent with the mesh PMLs."<<std::endl;
        std::cout<<"The number of PML types is "<<pmlTypeVec.size()<<" and there are "<<pmlMatIdVec.size();
        std::cout<<" PMLs in the mesh."<<std::endl;
        DebugStop();
    }
    if(boundMatIdVec.size() != boundTypeVec.size()){
        std::cout<<"The number of Boundary Conditions types is not consistent with the mesh."<<std::endl;
        std::cout<<"The number of BC types is "<<boundTypeVec.size()<<" and there are "<<boundMatIdVec.size();
        std::cout<<" boundaries in the mesh."<<std::endl;
        DebugStop();
    }
    const int outerMaterialPos = volMatIdVec.size() - 1;
    TPZCompMesh *cmeshH1 = new TPZCompMesh(gmesh);
    cmeshH1->SetDefaultOrder(pOrder); // seta ordem polimonial de aproximacao
    cmeshH1->SetDimModel(dim);        // seta dimensao do modelo
    // Inserindo material na malha
    const int nState = 1;
    TPZVec<STATE> sol; // only for creating material. this material will not be
    // used in reality

    TPZL2Projection *matH1;
    for (int i = 0; i < volMatIdVec.size(); ++i) {
        matH1 = new TPZL2Projection(volMatIdVec[i], dim, nState, sol);
        cmeshH1->InsertMaterialObject(matH1);
    }

    for (int i = 0; i < pmlMatIdVec.size(); ++i) {
        matH1 = new TPZL2Projection(pmlMatIdVec[i], dim, nState, sol);
        cmeshH1->InsertMaterialObject(matH1);
    }

    /// electrical and magnetic conductor boundary conditions
    TPZFNMatrix<1, STATE> val1(1, 1, 0.), val2(1, 1, 0.);
    TPZBndCond *bcond = nullptr;
    for(int i = 0; i < boundMatIdVec.size(); i++){
        bcond = matH1->CreateBC(matH1, boundMatIdVec[i], boundTypeVec[i], val1, val2);
        cmeshH1->InsertMaterialObject(bcond);
    }
    cmeshH1->AutoBuild();
    cmeshH1->CleanUpUnconnectedNodes();

    TPZCompMesh *cmeshHCurl = new TPZCompMesh(gmesh);
    cmeshHCurl->SetDefaultOrder(pOrder); // seta ordem polimonial de aproximacao
    cmeshHCurl->SetDimModel(dim);        // seta dimensao do modelo
    // Inserindo material na malha

    TPZVecL2 *matHCurl;
    for (int i = 0; i < volMatIdVec.size(); ++i) {
        matHCurl = new TPZVecL2(volMatIdVec[i]);
        cmeshHCurl->InsertMaterialObject(matHCurl);
    }

    for (int i = 0; i < pmlMatIdVec.size(); ++i) {
        matHCurl = new TPZVecL2(pmlMatIdVec[i]);
        cmeshHCurl->InsertMaterialObject(matHCurl);
    }

    for(int i = 0; i < boundMatIdVec.size(); i++){
        bcond = matHCurl->CreateBC(matHCurl, boundMatIdVec[i], boundTypeVec[i], val1, val2);
        cmeshHCurl->InsertMaterialObject(bcond);
    }

    cmeshHCurl->SetAllCreateFunctionsHCurl(); // define espaco de aproximacao

    cmeshHCurl->AutoBuild();
    cmeshHCurl->CleanUpUnconnectedNodes();

    TPZCompMesh *cmeshMF = new TPZCompMesh(gmesh);
    TPZManVector<TPZMatModalAnalysis *,8> matMultiPhysics( volMatIdVec.size() + pmlMatIdVec.size() );

    if (isCutOff) {
        for (int i = 0; i < volMatIdVec.size(); ++i) {
            matMultiPhysics[i] = new TPZMatWaveguideCutOffAnalysis(volMatIdVec[i], f0, urVec[i], erVec[i], 1./scale);
            cmeshMF->InsertMaterialObject(matMultiPhysics[i]);
        }
    } else {
        for (int i = 0; i < volMatIdVec.size(); ++i) {
            matMultiPhysics[i] = new TPZMatModalAnalysis(volMatIdVec[i], f0, urVec[i], erVec[i], 1./scale);
            cmeshMF->InsertMaterialObject(matMultiPhysics[i]);
        }
    }
    for(int i = 0; i < pmlMatIdVec.size(); i++){
        REAL xMax =-1e20,xMin = 1e20,yMax =-1e20,yMin =1e20;
        TPZGeoMesh * gmesh = cmeshMF->Reference();
        for (int iel = 0; iel < gmesh->NElements(); ++iel) {
            TPZGeoEl *geo = gmesh->Element(iel);
            if(geo->MaterialId() == pmlMatIdVec[i]){
                for (int iNode = 0; iNode < geo->NCornerNodes(); ++iNode) {
                    TPZManVector<REAL,3> co(3);
                    geo->Node(iNode).GetCoordinates(co);
                    const REAL & xP = co[0];
                    const REAL & yP = co[1];
                    if( xP > xMax ){
                        xMax = xP;
                    }
                    if( xP < xMin ){
                        xMin = xP;
                    }
                    if( yP > yMax ){
                        yMax = yP;
                    }
                    if( yP < yMin ){
                        yMin = yP;
                    }
                }
            }
        }
        bool attx, atty;
        REAL xBegin, yBegin, d;
        switch(pmlTypeVec[i]){
            case SPZModalAnalysisData::xp:
                attx = true; atty = false;
                xBegin = xMin; yBegin = -1;
                d = xMax - xMin;
                break;
            case SPZModalAnalysisData::yp:
                attx = false; atty = true;
                xBegin = -1; yBegin = yMin;
                d = yMax - yMin;
                break;
            case SPZModalAnalysisData::xm:
                attx = true; atty = false;
                xBegin = xMax; yBegin = -1;
                d = xMax - xMin;
                break;
            case SPZModalAnalysisData::ym:
                attx = false; atty = true;
                xBegin = -1; yBegin = yMax;
                d = yMax - yMin;
                break;
            case SPZModalAnalysisData::xpyp:
                attx = true; atty = true;
                xBegin = xMin; yBegin = yMin;
                d = xMax - xMin;
                break;
            case SPZModalAnalysisData::xmyp:
                attx = true; atty = true;
                xBegin = xMax; yBegin = yMin;
                d = xMax - xMin;
                break;
            case SPZModalAnalysisData::xmym:
                attx = true; atty = true;
                xBegin = xMax; yBegin = yMax;
                d = xMax - xMin;
                break;
            case SPZModalAnalysisData::xpym:
                attx = true; atty = true;
                xBegin = xMin; yBegin = yMax;
                d = xMax - xMin;
                break;
        }

        matMultiPhysics[volMatIdVec.size() + i] =
                new TPZMatWaveguidePml(pmlMatIdVec[i],
                                       *matMultiPhysics[outerMaterialPos],
                                       attx,xBegin,
                                       atty,yBegin,
                                       alphaMax,d);
        cmeshMF->InsertMaterialObject(matMultiPhysics[volMatIdVec.size() + i]);
    }


    for(int i = 0; i < boundMatIdVec.size(); i++){
        bcond = matMultiPhysics[0]->CreateBC(matMultiPhysics[0], boundMatIdVec[i], boundTypeVec[i], val1, val2);
        cmeshMF->InsertMaterialObject(bcond);
    }
    std::set<int> set;
    for (int i = 0; i < volMatIdVec.size(); ++i) {
        set.insert(volMatIdVec[i]);
    }
    for (int i = 0; i < pmlMatIdVec.size(); ++i) {
        set.insert(pmlMatIdVec[i]);
    }
    for (int i = 0; i < boundMatIdVec.size(); ++i) {
        set.insert(boundMatIdVec[i]);
    }

    cmeshMF->SetDimModel(dim);
    cmeshMF->SetAllCreateFunctionsMultiphysicElem();

    cmeshMF->AutoBuild(set);
    cmeshMF->CleanUpUnconnectedNodes();

    TPZVec<TPZCompMesh *> meshVec(2);
    meshVec[matMultiPhysics[0]->H1Index()] = cmeshH1;
    meshVec[matMultiPhysics[0]->HCurlIndex()] = cmeshHCurl;

    TPZBuildMultiphysicsMesh::AddElements(meshVec, cmeshMF);
    TPZBuildMultiphysicsMesh::AddConnects(meshVec, cmeshMF);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshVec, cmeshMF);

    cmeshMF->ExpandSolution();
    cmeshMF->ComputeNodElCon();
    cmeshMF->CleanUpUnconnectedNodes();

    if(print){
        std::ofstream fileH1(prefix + "cmeshH1.txt");
        cmeshH1->Print(fileH1);
        std::ofstream fileHCurl(prefix + "cmeshHCurl.txt");
        cmeshHCurl->Print(fileHCurl);
        std::ofstream fileMF(prefix + "cmeshMFHCurl.txt");
        cmeshMF->Print(fileMF);
    }

    meshVecOut.resize(3);

    meshVecOut[0] = cmeshMF;
    meshVecOut[1 + matMultiPhysics[0]->H1Index()] = cmeshH1;
    meshVecOut[1 + matMultiPhysics[0]->HCurlIndex()] = cmeshHCurl;
    return;
}
