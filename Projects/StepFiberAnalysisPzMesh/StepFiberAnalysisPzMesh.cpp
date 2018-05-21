/**
 * @file
 * @brief Performs modal analysis in an step-index optical fiber.
 * @details This project is aimed at electromagnetic analysis
 * of step-index optical fibers generating the geometric mesh in NeoPZ.
 * @author Francisco Orlandini
 * @since 2018
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
#include "tpzgeoblend.h"
#ifdef USING_SLEPC
#include <TPZSlepcEPSHandler.h>
#include <TPZSlepcSTHandler.h>
#include <Complex/TPZMatWaveguidePml.h>

#endif
#include "parameter_handler.h"
#include "TPZMatWaveguideCutOffAnalysis.h"
#include "TPZMatModalAnalysis.h"
#include "pzintel.h"
#include "SPZModalAnalysisDataReader.h"
#include "SPZModalAnalysisData.h"

void RunSimulation(SPZModalAnalysisData &simData,std::ostringstream & eigeninf, const int &hStep);

void CreateGMesh(TPZGeoMesh *&gmesh, const int &hStep, TPZVec<int> &matIdVec,
                 const bool &print, const REAL &scale = 1.);

void CreateCMesh(TPZVec<TPZCompMesh *> &meshVecOut, TPZGeoMesh *gmesh,
                 int pOrder, TPZVec<int> &matIdVec, TPZVec<STATE> &urVec,
                 TPZVec<STATE> &erVec, REAL f0, bool isCutOff, const std::string &prefix,
                 const bool &print, const REAL &scale, const bool &hasPML, const REAL &alphaMax);

void FilterBoundaryEquations(TPZVec<TPZCompMesh *> cmeshMF,
                             TPZVec<long> &activeEquations, int &neq,
                             int &neqOriginal);


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
            simData.pzOpts.pOrder = pOrderOrig;
            for (int iP = 0; iP < simData.pzOpts.pSteps; ++iP) {
                std::cout<<"freq step: "<<iFreq+1<<" out of "<<simData.physicalOpts.freqVec.size()<<std::endl;
                std::cout<<"h    step: "<< iH+1<<" out of "<<simData.pzOpts.hSteps<<std::endl;
                std::cout<<"Beginning p step "<<iP+1<<" out of "<<simData.pzOpts.pSteps<<std::endl;
                std::ostringstream eigeninfo;
                RunSimulation(simData,eigeninfo, iH);
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

void RunSimulation(SPZModalAnalysisData &simData,std::ostringstream &eigeninfo, const int &hStep) {
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    std::cout<<"Creating GMesh..."<<std::endl;
    boost::posix_time::ptime t1_g =
        boost::posix_time::microsec_clock::local_time();
    TPZVec<int> matIdVec;
    CreateGMesh(gmesh, hStep, matIdVec, simData.pzOpts.exportGMesh, simData.pzOpts.scaleFactor);
    boost::posix_time::ptime t2_g =
        boost::posix_time::microsec_clock::local_time();
    std::cout<<"Created!  "<<t2_g-t1_g<<std::endl;
    TPZVec<TPZCompMesh *> meshVec(1);
    std::cout<<"Creating CMesh..."<<std::endl;
    boost::posix_time::ptime t1_c =
        boost::posix_time::microsec_clock::local_time();

    CreateCMesh(meshVec, gmesh, simData.pzOpts.pOrder, matIdVec,
                simData.physicalOpts.urVec, simData.physicalOpts.erVec,
                simData.physicalOpts.lambda,simData.physicalOpts.isCutOff,
                simData.pzOpts.prefix,simData.pzOpts.exportCMesh,
                simData.pzOpts.scaleFactor,
                simData.physicalOpts.hasPML, simData.physicalOpts.alphaMax); // funcao para criar a malha computacional
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
        an.PostProcess(3);
    }
//    std::cout << "Assembling..." << std::endl;
//    boost::posix_time::ptime t1 =
//        boost::posix_time::microsec_clock::local_time();
//    an.Assemble();
//    boost::posix_time::ptime t2 =
//        boost::posix_time::microsec_clock::local_time();
//    std::cout << "Finished assembly." << std::endl;

//     std::cout << "Solving..." << std::endl;
//     boost::posix_time::ptime t3 =
//         boost::posix_time::microsec_clock::local_time();

//     an.Solve();
//     boost::posix_time::ptime t4 =
//         boost::posix_time::microsec_clock::local_time();
//     std::cout << "Time for assembly " << t2 - t1 << " Time for solving "
//               << t4 - t3 << std::endl;
//     const TPZManVector<SPZAlwaysComplex<STATE>::type> eigenValues = an.GetEigenvalues();
//     const TPZFMatrix<SPZAlwaysComplex<STATE>::type> eigenVectors = an.GetEigenvectors();

//     typedef std::numeric_limits< double > dbl;
//     std::cout.precision(dbl::max_digits10);
//     for(int i = 0; i < eigenValues.size() ; i++){
//         std::cout<<std::fixed<<eigenValues[i]<<std::endl;
//     }
//     if(simData.pzOpts.exportEigen){
//         eigeninfo.precision(dbl::max_digits10);
//         std::cout<<"Exporting eigen info..."<<std::endl;
//         REAL hSize = -1e12;
//         REAL elRadius = 0;
//         TPZVec<REAL> qsi(2,0.25);
//         TPZVec<REAL> x(3,0.);
//         for (int j = 0; j < gmesh->NElements(); ++j) {
//             TPZGeoEl &el = *(gmesh->ElementVec()[j]);
//             // for (int i = 0; i < el.NCornerNodes() ; ++i) {
//             //     auto node = el.Node(i);
//             //     if(node.Coord(0) < tol && node.Coord(1) < tol){
//             //         elRadius = el.ElementRadius();
//             //         hSize = elRadius < hSize ? elRadius : hSize;
//             //     }
//             // }
//            TPZBndCond *matBound = dynamic_cast<TPZBndCond *>(meshVec[0]->MaterialVec()[el.MaterialId()]);
//            TPZMatWaveguidePml *matPml = dynamic_cast<TPZMatWaveguidePml *>(meshVec[0]->MaterialVec()[el.MaterialId()]);
//            if(!matBound && !matPml){
//                elRadius = el.ElementRadius();
//                hSize = elRadius > hSize ? elRadius : hSize;
//            }
//         }
//         eigeninfo << std::fixed << neq << "," << gmesh->NElements() << ",";
//         eigeninfo << std::fixed << hSize << "," << simData.pzOpts.pOrder<<",";
//         eigeninfo << std::fixed << simData.physicalOpts.lambda<<",";
//         eigeninfo << eigenValues.size()<<",";
//         for(int i = 0; i < eigenValues.size() ; i++){
//             eigeninfo<<std::fixed<<std::real(eigenValues[i])<<",";
//             eigeninfo<<std::fixed<<std::imag(eigenValues[i]);
//             if(i != eigenValues.size() - 1 ) {
//                 eigeninfo << ",";
//             }
//         }
//         eigeninfo << std::endl;

//         std::cout<<" Done!"<<std::endl;
//     }

//     if (simData.pzOpts.genVTK) {
//         TPZMatModalAnalysis *matPointer =
//                 dynamic_cast<TPZMatModalAnalysis *>(meshVec[0]->MaterialVec()[1]);
//         TPZVec<TPZCompMesh *> temporalMeshVec(2);
//         temporalMeshVec[matPointer->H1Index()] = meshVec[1 + matPointer->H1Index()];
//         temporalMeshVec[matPointer->HCurlIndex()] =
//                 meshVec[1 + matPointer->HCurlIndex()];

//         std::cout << "Post Processing..." << std::endl;

//         TPZStack<std::string> scalnames, vecnames;
//         scalnames.Push("Ez");
//         vecnames.Push("Et");
//         std::string plotfile = simData.pzOpts.prefix + "fieldPlot" + ".vtk";
//                                                         // estara na pasta debug
//         const int dim = 2;
//         an.DefineGraphMesh(dim, scalnames, vecnames,
//                            plotfile);  // define malha grafica

//         TPZFMatrix<SPZAlwaysComplex<STATE>::type> currentEigenvector(neq,1);
//         TPZFMatrix<SPZAlwaysComplex<STATE>::type> scatteredEigen(neqOriginal,1);
//         for (int iSol = 0; iSol < eigenValues.size(); iSol++) {
//             for(int j = 0; j < eigenVectors.Rows(); j++){
//                 currentEigenvector(j,0) = simData.pzOpts.absVal ?
//                                           std::abs(eigenVectors.GetVal(j,iSol)) :
//                                           std::real(eigenVectors.GetVal(j,iSol));
//             }
//             strmtrx->EquationFilter().Scatter(currentEigenvector, scatteredEigen);
//             an.LoadSolution(scatteredEigen);
//             TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(temporalMeshVec,
//                                                                cmesh);
//             an.PostProcess(simData.pzOpts.vtkRes);
//         }
//     }
//     gmesh->SetReference(nullptr);
//     for (int k = 0; k < meshVec.size(); ++k) {
//         meshVec[k]->SetReference(nullptr);
//         delete meshVec[k];
//         meshVec[k] = nullptr;
//     }
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



void
CreateGMesh(TPZGeoMesh *&gmesh, const int &hStep, TPZVec<int> &matIdVec,
          const bool &print, const REAL &scale) {
    const int nDivTCore = 3, nDivRCore = 3 , nDivTCladding = 3;

    const int matIdCore = 1, matIdCladding = 2, matIdBC= -1;
    const int matIdPMLxp = 91,
              matIdPMLyp = 92,
              matIdPMLxm = 93,
              matIdPMLym = 94,
              matIdPMLxpyp = 95,
              matIdPMLxmyp = 96,
              matIdPMLxmym = 97,
              matIdPMLxpym = 98;

    const REAL rCore = 1.;
    std::vector<REAL> xc(2);
    xc[0] = 0.;
    xc[1] = 0.;
    const REAL bound = 2.;
    TPZVec<REAL> ptCircle_1(2,0.);
    ptCircle_1[0] = xc[0] + rCore * cos(1*M_PI/4);
    ptCircle_1[1] = xc[1] + rCore * sin(1*M_PI/4);
    TPZVec<REAL> ptCircle_2(2,0.);
    ptCircle_2[0] = xc[0] + rCore * cos(3*M_PI/4);
    ptCircle_2[1] = xc[1] + rCore * sin(3*M_PI/4);
    TPZVec<REAL> ptCircle_3(2,0.);
    ptCircle_3[0] = xc[0] + rCore * cos(5*M_PI/4);
    ptCircle_3[1] = xc[1] + rCore * sin(5*M_PI/4);
    TPZVec<REAL> ptCircle_4(2,0.);
    ptCircle_4[0] = xc[0] + rCore * cos(7*M_PI/4);
    ptCircle_4[1] = xc[1] + rCore * sin(7*M_PI/4);

    TPZVec<REAL> ptSquare_1(2,0.);
    ptSquare_1[0] = ptCircle_1[0] / 2.;
    ptSquare_1[1] = ptCircle_1[1] / 2.;
    TPZVec<REAL> ptSquare_2(2,0.);
    ptSquare_2[0] = ptCircle_2[0] / 2.;
    ptSquare_2[1] = ptCircle_2[1] / 2.;
    TPZVec<REAL> ptSquare_3(2,0.);
    ptSquare_3[0] = ptCircle_3[0] / 2.;
    ptSquare_3[1] = ptCircle_3[1] / 2.;
    TPZVec<REAL> ptSquare_4(2,0.);
    ptSquare_4[0] = ptCircle_4[0] / 2.;
    ptSquare_4[1] = ptCircle_4[1] / 2.;

    TPZVec<REAL> ptBound_1(2,0.);
    ptBound_1[0] = 1 * bound;
    ptBound_1[1] = 1 * bound;
    TPZVec<REAL> ptBound_2(2,0.);
    ptBound_2[0] = -1 * bound;
    ptBound_2[1] = 1 * bound;
    TPZVec<REAL> ptBound_3(2,0.);
    ptBound_3[0] = -1 * bound;
    ptBound_3[1] = -1 * bound;
    TPZVec<REAL> ptBound_4(2,0.);
    ptBound_4[0] = 1 * bound;
    ptBound_4[1] = -1 * bound;

    auto map_quad_side_arc = [](const TPZVec<REAL> &theta ,const TPZVec<REAL> &xc, const REAL &r, const REAL & s) {
        TPZVec<REAL> point(2,0.);
        point[0] = xc[0] + r*cos((theta[1] + s*theta[1] + theta[0] - s*theta[0])/2.);
        point[1] = xc[1] + r*sin((theta[1] + s*theta[1] + theta[0] - s*theta[0])/2.);
        return point;
    };


    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    // quadrilateral one - center square                                                                    //
    // quadrilateral two to five - fiber's core, from the right, counter-clockwise                          //
    // quadrilateral six to nine - fiber's cladding, from the right, counter-clockwise                      //
    // quadrilateral ten to thirteen - PML regions (excluding corners), from the right, counter-clockwise   //
    // quadrilateral fourteen to seventeen - PML corners, from the top-right, counter-clockwise             //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////

    typedef struct QuadrilateralData {
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

        TPZFMatrix<REAL> side1pts;
        TPZFMatrix<REAL> side2pts;
        TPZFMatrix<REAL> side3pts;
        TPZFMatrix<REAL> side4pts;
        TPZFMatrix<REAL> interiorPts;

        QuadrilateralData(const TPZVec<REAL> &c1,const TPZVec<REAL> &c2,const TPZVec<REAL> &c3,const TPZVec<REAL> &c4) :
                coord1(c1),coord2(c2),coord3(c3),coord4(c4){
            auto map_quad_side_linear = [](const TPZVec<REAL> &coord1 ,const TPZVec<REAL> &coord2,const REAL & s) {
                TPZVec<REAL> point(2,0.);
                point[0] = (coord1[0] - s*coord1[0] + coord2[0] + s*coord2[0])/2.;
                point[1] = (coord1[1] - s*coord1[1] + coord2[1] + s*coord2[1])/2.;
                return point;
            };
            mapSide1 = std::bind(map_quad_side_linear,coord1,coord2,std::placeholders::_1);
            mapSide2 = std::bind(map_quad_side_linear,coord2,coord3,std::placeholders::_1);
            mapSide3 = std::bind(map_quad_side_linear,coord3,coord4,std::placeholders::_1);
            mapSide4 = std::bind(map_quad_side_linear,coord4,coord1,std::placeholders::_1);

        }
        void SetMapsLinearity(const bool &m1, const bool &m2, const bool &m3, const bool &m4){
            isSide1nonLinear = m1;
            isSide2nonLinear = m2;
            isSide3nonLinear = m3;
            isSide4nonLinear = m4;
        }

        void CreateQuadrilateral(const int &nQsi, const int &nEta){

            const int nPtsSide1 = nQsi + (nQsi - 1)*isSide1nonLinear;
            const int nPtsSide2 = nEta + (nEta - 1)*isSide2nonLinear;
            const int nPtsSide3 = nQsi + (nQsi - 1)*isSide3nonLinear;
            const int nPtsSide4 = nEta + (nEta - 1)*isSide4nonLinear;

            const int nPtsInterior = (nEta - 2) * (nQsi - 2);

            // std::cout<<"nPtsSide1 :"<<nPtsSide1<<" ";
            // std::cout<<"nPtsSide2 :"<<nPtsSide2<<" ";
            // std::cout<<"nPtsSide3 :"<<nPtsSide3<<" ";
            // std::cout<<"nPtsSide4 :"<<nPtsSide4<<" ";
            // std::cout<<"nPtsInterior :"<<nPtsInterior<<" "<<std::endl;

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
                    x[0] = coord1side[1] * (1.-s)/2 + coord2side[1] * (1+s)/2;
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

                        std::cout<<"x:   "<<sidePts(iPt,0)<<", y:   "<<sidePts(iPt,1)<<std::endl;
                        iPt++;
                    }
                }
            };

            std::cout<<"---------------side 1---------------"<<std::endl;
            getPoints(side1pts,nPtsSide1,1,-1., 1.,-1.,-1.);
            std::cout<<"---------------side 2---------------"<<std::endl;
            getPoints(side2pts,1,nPtsSide2, 1., 1.,-1., 1.);
            std::cout<<"---------------side 3---------------"<<std::endl;
            getPoints(side3pts,nPtsSide3,1, 1.,-1., 1., 1.);
            std::cout<<"---------------side 4---------------"<<std::endl;
            getPoints(side4pts,1,nPtsSide4,-1.,-1., 1.,-1.);
            const REAL deltaQsiInt = 2/(nQsi - 1);
            const REAL deltaEtaInt = 2/(nEta - 1);
            //std::cout<<"--------------- intP ---------------"<<std::endl;
            getPoints(interiorPts,nQsi-2,nEta-2,-1. + deltaQsiInt, 1. - deltaQsiInt,-1. + deltaEtaInt, 1. - deltaEtaInt );
        }
    } QuadrilateralData;

    TPZVec<QuadrilateralData *> quadVec(17,nullptr);

    quadVec[0] = new QuadrilateralData(ptSquare_1,ptSquare_2,ptSquare_3,ptSquare_4);
    quadVec[0]->SetMapsLinearity(false,false,false,false);
    quadVec[0]->CreateQuadrilateral(nDivRCore,nDivRCore);

    TPZVec<REAL> theta(2,0.);
    theta[0] =-1 * M_PI / 4;
    theta[1] =1 * M_PI / 4;

    quadVec[1] = new QuadrilateralData(ptCircle_1,ptSquare_1,ptSquare_3,ptCircle_4);
    quadVec[1]->mapSide4 = std::bind(map_quad_side_arc,theta, xc, rCore,std::placeholders::_1);
    quadVec[1]->SetMapsLinearity(false,false,false,true);
    quadVec[1]->CreateQuadrilateral(nDivTCore,nDivRCore);
    //std::function<std::vector<REAL>(const REAL &)> el_1_side_1 = std::bind(map_quad_side_linear,ptBound_3,ptCircle_3,std::placeholders::_1);
    //std::function<std::vector<REAL>(const REAL &)> el_1_side_2 = std::bind(map_quad_side_arc,theta,xc,rCore,std::placeholders::_1);
    //std::function<std::vector<REAL>(const REAL &)> el_1_side_3 = std::bind(map_quad_side_linear,ptCircle_4,ptBound_2,std::placeholders::_1);
    //std::function<std::vector<REAL>(const REAL &)> el_1_side_4 = std::bind(map_quad_side_linear,ptBound_2,ptBound_3,std::placeholders::_1);





    //CreateQuadrilateral(ptBound_3,ptCircle_3,ptCircle_2,ptBound_2,el_1_side_1,el_1_side_2,el_1_side_3,el_1_side_4,false,true,false,false,nx,ny);



    ///etc etc
    matIdVec.Resize(3);
    matIdVec[0]=1;
    matIdVec[1]=2;
    matIdVec[2]=-1;
    const int nNodes = 5;

    gmesh = new TPZGeoMesh;
    gmesh->NodeVec().Resize(nNodes);
    //creates core
    {

        // create 4 nodes which will be triangle vertices
        // r arg 0, r arg pi/2, r arg pi, r arg 3pi/2
        int nodeId = 0;
        for (int iNode = 0; iNode < 4; iNode++) {
            TPZManVector<REAL, 3> node(3, 0.);
            const int c0 = -1 + (int)2*(int)(((iNode+1)%4)/2);//expected: -1 1 1 -1
            const int c1 = -1 + 2* (iNode / 2);//expected: -1 -1 1 1
            node[0] = c0*rCore;
            node[1] = c1*rCore;
            gmesh->NodeVec()[nodeId].Initialize(node, *gmesh);
            nodeId++;
        }

        TPZManVector<REAL, 3> node(3, 0.);
        node[0] = 0.5;
        node[1] = -0.5;
        gmesh->NodeVec()[nodeId].Initialize(node, *gmesh);

        TPZManVector<long, 3> nodesIdVec(3, 0);
        // creates volumetric elements
        TPZManVector<int,2> initialNodes(2,0);
        initialNodes[0] = 0;
        initialNodes[1] = 2;

        for (int iTri = 0; iTri < 2; iTri++) {
            nodesIdVec[0] = initialNodes[iTri];
            nodesIdVec[1] = (initialNodes[iTri] + 1) % 4;
            nodesIdVec[2] = (initialNodes[iTri] + 2) % 4;

            TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle> >  * triangulo =
                  new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle> >  (nodesIdVec,iTri+1,*gmesh);
        }
        nodesIdVec[0] = 0;
        nodesIdVec[1] = 2;
        nodesIdVec[2] = 4;//midpoint
        TPZGeoElRefPattern< pzgeom::TPZArc3D > *arc =
                new TPZGeoElRefPattern< pzgeom::TPZArc3D > (nodesIdVec,2,*gmesh);
    }

    gmesh->BuildConnectivity();
    if(print){
        std::string meshFileName = "gmesh";
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


void CreateCMesh(TPZVec<TPZCompMesh *> &meshVecOut, TPZGeoMesh *gmesh,
                 int pOrder, TPZVec<int> &matIdVec, TPZVec<STATE> &urVec,
                 TPZVec<STATE> &erVec, REAL f0, bool isCutOff, const std::string &prefix,
                 const bool &print, const REAL &scale, const bool &hasPML, const REAL &alphaMax) {
    enum {
      dirichlet = 0
    }; // tipo da condicao de contorno do problema
    const int dim = 2;   // dimensao do problema

    TPZVec<int> volMatId(matIdVec.size()-1);
    for(int i = 0; i< volMatId.size(); i++){
      volMatId[i] = matIdVec[i];
    }
    const int boundMatId = matIdVec[matIdVec.size()-1];
    const int nPML = hasPML ? 8 : 0;
    if(volMatId.size() - nPML !=urVec.size()) DebugStop();

    /// criar malha computacional H1
    TPZCompMesh *cmeshH1 = new TPZCompMesh(gmesh);
    cmeshH1->SetDefaultOrder(pOrder); // seta ordem polimonial de aproximacao
    cmeshH1->SetDimModel(dim);        // seta dimensao do modelo
    // Inserindo material na malha
    const int nState = 1;
    TPZVec<STATE> sol; // only for creating material. this material will not be
    // used in reality

    TPZL2Projection *matH1;
    for (int i = 0; i < volMatId.size(); ++i) {
        matH1 = new TPZL2Projection(volMatId[i], dim, nState, sol);
        cmeshH1->InsertMaterialObject(matH1);
    }

    /// electrical conductor boundary conditions
    TPZFNMatrix<1, STATE> val1(1, 1, 0.), val2(1, 1, 0.);

    TPZMaterial *BCondH1Dir = matH1->CreateBC(
        matH1, boundMatId, dirichlet, val1, val2); // cria material que implementa a
    // condicao de contorno de dirichlet

    cmeshH1->InsertMaterialObject(BCondH1Dir); // insere material na malha
    // Cria elementos computacionais que gerenciarao o espaco de aproximacao da
    // malha
    cmeshH1->AutoBuild();
    cmeshH1->CleanUpUnconnectedNodes();

    TPZCompMesh *cmeshHCurl = new TPZCompMesh(gmesh);
    cmeshHCurl->SetDefaultOrder(pOrder); // seta ordem polimonial de aproximacao
    cmeshHCurl->SetDimModel(dim);        // seta dimensao do modelo
    // Inserindo material na malha

    TPZVecL2 *matHCurl;
    for (int i = 0; i < volMatId.size(); ++i) {
        matHCurl = new TPZVecL2(volMatId[i]);
        cmeshHCurl->InsertMaterialObject(matHCurl);
    }

    val1(0, 0) = 0.;
    val2(0, 0) = 0.;
    TPZMaterial *BCondHCurlDir =
        matHCurl->CreateBC(matHCurl, boundMatId, dirichlet, val1,
                           val2); // cria material que implementa a condicao de
    // contorno de dirichlet

    cmeshHCurl->InsertMaterialObject(BCondHCurlDir); // insere material na malha

    cmeshHCurl->SetAllCreateFunctionsHCurl(); // define espaco de aproximacao

    cmeshHCurl->AutoBuild();
    cmeshHCurl->CleanUpUnconnectedNodes();

    TPZCompMesh *cmeshMF = new TPZCompMesh(gmesh);
    TPZVec<TPZMatModalAnalysis *>matMultiPhysics(volMatId.size());

    ///PML Materials
    TPZVec<int> pmlIds(0,0);
    int outerMaterial = 0;
    if(hasPML){
        pmlIds.Resize(8);
        for(int i = 0; i < 8 ; i++){
            pmlIds[i] = volMatId[volMatId.size()-8 + i];
        }
        volMatId.Resize(volMatId.size()-8);
        outerMaterial = volMatId.size()-1;

        // for (int i = 0; i < cmeshH1->NElements(); ++i) {
        //     TPZCompEl *compel = cmeshH1->Element(i);
        //     TPZGeoEl *geo = cmeshH1->Element(i)->Reference();
        //     const int matId = geo->MaterialId();
        //     if (matId == pmlIds[0] ||
        //         matId == pmlIds[1] ||
        //         matId == pmlIds[2] ||
        //         matId == pmlIds[3] ||
        //         matId == pmlIds[4] ||
        //         matId == pmlIds[5] ||
        //         matId == pmlIds[6] ||
        //         matId == pmlIds[7]) {
        //         TPZInterpolatedElement *cel = dynamic_cast<TPZInterpolatedElement *>(compel);
        //         cel->PRefine(2);
        //     }
        // }
        // for (int i = 0; i < cmeshHCurl->NElements(); ++i) {
        //     TPZCompEl *compel = cmeshHCurl->Element(i);
        //     TPZGeoEl *geo = cmeshHCurl->Element(i)->Reference();
        //     const int matId = geo->MaterialId();
        //     if (matId == pmlIds[0] ||
        //         matId == pmlIds[1] ||
        //         matId == pmlIds[2] ||
        //         matId == pmlIds[3] ||
        //         matId == pmlIds[4] ||
        //         matId == pmlIds[5] ||
        //         matId == pmlIds[6] ||
        //         matId == pmlIds[7]) {
        //         TPZInterpolatedElement *cel = dynamic_cast<TPZInterpolatedElement *>(compel);
        //         cel->PRefine(2);
        //     }
        // }

        // cmeshH1->ExpandSolution();
        // cmeshH1->ComputeNodElCon();
        // cmeshH1->CleanUpUnconnectedNodes();

        // cmeshHCurl->ExpandSolution();
        // cmeshHCurl->ComputeNodElCon();
        // cmeshHCurl->CleanUpUnconnectedNodes();

    }

    if (isCutOff) {
        for (int i = 0; i < volMatId.size(); ++i) {
            matMultiPhysics[i] =
                new TPZMatWaveguideCutOffAnalysis(volMatId[i], f0, urVec[i], erVec[i], 1./scale);
            cmeshMF->InsertMaterialObject(matMultiPhysics[i]);
        }
    } else {
        for (int i = 0; i < volMatId.size(); ++i) {
            matMultiPhysics[i] =
                new TPZMatModalAnalysis(volMatId[i], f0, urVec[i], erVec[i], 1./scale);
            cmeshMF->InsertMaterialObject(matMultiPhysics[i]);
        }
    }

    if(hasPML){
        REAL xR=1e20,xL=-1e20,yU=1e20,yD=-1e20,d=-1e20;
        TPZGeoMesh * gmesh = cmeshMF->Reference();
        const int urId = pmlIds[5];
        const int llId = pmlIds[7];
        for (int i = 0; i < gmesh->NElements(); ++i) {
            TPZGeoEl *geo = gmesh->Element(i);
            if(geo->MaterialId() == urId){
                for (int j = 0; j < geo->NCornerNodes(); ++j) {
                    TPZManVector<REAL,3> co(3);
                    geo->Node(j).GetCoordinates(co);
                    const REAL & xP = co[0];
                    const REAL & yP = co[1];
                    if( xP < xR ){
                        xR = xP;
                    }
                    if( yP < yU ){
                        yU = yP;
                    }
                    if( xP > d){
                        d = xP;
                    }
                }
            }
            else if(geo->MaterialId() == llId){
                for (int j = 0; j < geo->NCornerNodes(); ++j) {
                    TPZManVector<REAL,3> co(3);
                    geo->Node(j).GetCoordinates(co);
                    const REAL & xP = co[0];
                    const REAL & yP = co[1];
                    if( xP > xL ){
                        xL = xP;
                    }
                    if( yP > yD ){
                        yD = yP;
                    }
                }
            }
            else{
                continue;
            }
        }
        d = d - xR;
        REAL xP;
        REAL yP;
        bool attx = true;
        bool atty = false;
        for(int i = 0; i < nPML; i++){
            xP = (i / 2) % 2 ? xL : xR;
            yP = !((i / 2) % 2) != !(i%2)? yU : yD;
            if(i <4){
                attx = !(i % 2);
                atty = !(attx);
            }
            if(i == 4){
                attx = true;
                atty = true;
            }
            matMultiPhysics[volMatId.size() + i] =
                    new TPZMatWaveguidePml(pmlIds[i],
                                           *matMultiPhysics[outerMaterial],
                                           attx,xP,
                                           atty,yP,
                                           alphaMax,d);
            cmeshMF->InsertMaterialObject(matMultiPhysics[volMatId.size() + i]);
        }
    }

    val1(0, 0) = 0.;
    val2(0, 0) = 0.;
    TPZMaterial *BCondMFDir =
        matMultiPhysics[0]->CreateBC(matMultiPhysics[0], boundMatId, dirichlet, val1,
                                     val2); // any material is fine

    cmeshMF->InsertMaterialObject(BCondMFDir); // insere material na malha
    std::set<int> set;
    for (int i = 0; i < volMatId.size(); ++i) {
        set.insert(volMatId[i]);
    }
    for (int i = 0; i < nPML; ++i) {
        set.insert(pmlIds[i]);
    }
    set.insert(boundMatId);

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
