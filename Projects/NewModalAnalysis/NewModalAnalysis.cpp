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
#include <tpzgeoblend.h>
#include "pzgengrid.h"
#include "pzgmesh.h"
#include "TPZVTKGeoMesh.h"
#include "pzbndcond.h"
#include "pzlog.h"
#include "pzstrmatrix.h"
#include "pzl2projection.h"
#include "pzfstrmatrix.h"
#include "pzbuildmultiphysicsmesh.h"
#include "tpzgeoelrefpattern.h"
#include "boost/date_time/posix_time/posix_time.hpp"
#ifdef USING_SLEPC
#include <TPZSlepcEPSHandler.h>
#include <TPZSlepcSTHandler.h>
#include <Complex/TPZMatWaveguidePml.h>
#include <Complex/TPZMatWaveguidePmlHDiv.h>
#include <tpzgeoblend.h>

#endif
#include "parameter_handler.h"
#include "TPZMatWaveguideCutOffAnalysis.h"
#include "TPZMatModalAnalysis.h"
#include "TPZMatModalAnalysisHDiv.h"
#include "pzintel.h"
#include "TPZGmshReader.h"
#include "SPZModalAnalysisDataReader.h"
#include "SPZModalAnalysisData.h"

void RunSimulation(SPZModalAnalysisData &simData,std::ostringstream &eigeninfo, const int &factor);

typedef struct QuadrilateralData QuadrilateralData;
typedef struct EdgeData EdgeData;

TPZGeoMesh* CreateStructuredMesh(const TPZVec<TPZVec<REAL>> &pointsVec, TPZVec<EdgeData> &edgesVec, const TPZFMatrix<int> &quadPointsVec,
                                 const TPZVec<int> &matIdsQuads, const TPZVec<int> &nDivQsi, const TPZVec<int> &nDivEta,
                                 const TPZVec<bool> &side1NonLinearVec, const TPZVec<bool> &side2NonLinearVec,
                                 const TPZVec<bool> &side3NonLinearVec, const TPZVec<bool> &side4NonLinearVec,
                                 const TPZVec<TPZVec<REAL>> &thetaVec, const TPZVec<TPZVec<REAL> *> &xcRef,
                                 const TPZVec<int> &matIdBoundVec, const TPZVec<REAL> &boundDistVec,
                                 const TPZVec<REAL> &rVec, const bool nonLinearMapping = true);
void
CreateGMeshRectangularWaveguide(TPZGeoMesh *&gmesh, const std::string mshFileName, TPZVec<int> &matIdVec,
                                const std::string &prefix, const bool &print, const REAL &scale, const int &factor,
                                const REAL wDomain, const REAL hDomain, TPZVec<SPZModalAnalysisData::boundtype> &boundTypeVec,
                                TPZVec<SPZModalAnalysisData::pmltype> &pmlTypeVec, const bool &usingSymmetry,
                                const SPZModalAnalysisData::boundtype &symmetryType,
                                bool refine = false);
void
CreateStepFiberMesh(TPZGeoMesh *&gmesh, const std::string mshFileName, TPZVec<int> &matIdVec, const std::string &prefix,
                    const bool &print, const REAL &scale, const int &factor,
                    TPZVec<SPZModalAnalysisData::boundtype> &boundTypeVec,
                    TPZVec<SPZModalAnalysisData::pmltype> &pmlTypeVec, const REAL &realRCore, const REAL &dPML,
                    const REAL &boundDist, const REAL &outerReffIndex, const int &nLayersPml,
                    bool refine=false);

void
CreateHoleyFiberMesh(TPZGeoMesh *&gmesh, const std::string mshFileName, TPZVec<int> &matIdVec,
                     const std::string &prefix, const bool &print, const REAL &scale, const int &factor,
                     TPZVec<SPZModalAnalysisData::boundtype> &boundTypeVec,
                     TPZVec<SPZModalAnalysisData::pmltype> &pmlTypeVec, const REAL &dPML, const REAL &boundDist,
                     SPZModalAnalysisData::boundtype symmetryX, SPZModalAnalysisData::boundtype symmetryY,
                     bool refine=false,
                     TPZVec<std::function<bool (const TPZVec<REAL> &)>> refineRules =
                             TPZVec<std::function<bool (const TPZVec<REAL> &)>>(0,nullptr));

void
ReadGMesh(TPZGeoMesh *&gmesh, const std::string mshFileName, TPZVec<int> &matIdVec, const std::string &prefix, const bool &print, const REAL &scale = 1.);

void CreateCMesh(TPZVec<TPZCompMesh *> &meshVecOut, TPZGeoMesh *gmesh, int pOrder, const TPZVec<int> &matIdVec,
                 TPZVec<STATE> &urVec, TPZVec<STATE> &erVec, REAL f0, bool isCutOff, const std::string &prefix,
                 const bool &print, const REAL &scale, const REAL &alphaMax,
                 TPZVec<SPZModalAnalysisData::boundtype> &boundTypeVec,
                 TPZVec<SPZModalAnalysisData::pmltype> &pmlTypeVec, bool refineP = false,
                 TPZVec<std::function<bool (const TPZVec<REAL> &)>> refineRules =
                         TPZVec<std::function<bool (const TPZVec<REAL> &)>>(0,nullptr),
                         SPZModalAnalysisData::NedEl elType = SPZModalAnalysisData::NedEl::TypeOne);

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
                        simData.pzOpts.factorVec[iH] = iH;
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
    bool refineH = false;
    bool refineP = false;
    TPZVec<std::function<bool (const TPZVec<REAL> &)>> refineRulesH(0,nullptr);
    TPZVec<std::function<bool (const TPZVec<REAL> &)>> refineRulesP(0,nullptr);
    if(simData.pzOpts.usingNeoPzMesh){
        switch(simData.pzOpts.pzCase){

            case SPZModalAnalysisData::HoleyFiber:{
                refineH = simData.physicalOpts.holeyFiberOpts.refineH;
                refineP = simData.physicalOpts.holeyFiberOpts.refineP;
                const REAL dPML = simData.physicalOpts.holeyFiberOpts.dPML;
                const REAL boundDist = simData.physicalOpts.holeyFiberOpts.boundDist;
                const SPZModalAnalysisData::boundtype symmetryX = simData.physicalOpts.holeyFiberOpts.symmetryX;
                const SPZModalAnalysisData::boundtype symmetryY = simData.physicalOpts.holeyFiberOpts.symmetryY;

                const REAL rCore = simData.pzOpts.scaleFactor * 2.5e-6;
                const REAL Lambda = simData.pzOpts.scaleFactor * 6.75e-6;
                //holes centers
                TPZManVector<REAL,2> xc1(2, 0.);
                TPZManVector<REAL,2> xc2(2, 0.);
                xc1[0] = Lambda * 0.5;
                xc1[1] = Lambda * std::sqrt(3) * 0.5;
                xc2[0] = Lambda;
                xc2[1] = 0;
                REAL bDist = simData.physicalOpts.holeyFiberOpts.boundDist *simData.pzOpts.scaleFactor;
                if(refineP == true){
                    auto firstRuleHFiberP = [xc1,xc2,rCore,bDist,Lambda](const TPZVec<REAL> &xPos) {
                        return (
                                sqrt(xPos[0]*xPos[0]+xPos[1]*xPos[1])<(Lambda + rCore)
//                                &&
//                                sqrt((xPos[0] - xc1[0]) * (xPos[0] - xc1[0]) +
//                                     (xPos[1] - xc1[1]) * (xPos[1] - xc1[1])) > 0.6* rCore
//                                &&
//                                sqrt((xPos[0] - xc2[0]) * (xPos[0] - xc2[0]) +
//                                     (xPos[1] - xc2[1]) * (xPos[1] - xc2[1])) > 0.6* rCore
                        );
                    };
                    auto secondRuleHFiberP = [xc1,xc2,rCore,bDist](const TPZVec<REAL> &xPos) {
                        return (
                                sqrt(xPos[0]*xPos[0]+xPos[1]*xPos[1])<0.78
                                                                      *bDist
//                                &&
//                                sqrt((xPos[0] - xc1[0]) * (xPos[0] - xc1[0]) +
//                                     (xPos[1] - xc1[1]) * (xPos[1] - xc1[1])) > 0.6* rCore
//                                &&
//                                sqrt((xPos[0] - xc2[0]) * (xPos[0] - xc2[0]) +
//                                     (xPos[1] - xc2[1]) * (xPos[1] - xc2[1])) > 0.6* rCore
                        );
                    };
//                    auto thirdRuleHFiberP = [xc1,xc2,rCore,bDist](const TPZVec<REAL> &xPos) {
//                        return (
//                                sqrt(xPos[0]*xPos[0]+xPos[1]*xPos[1])<0.9*bDist
////                                &&
////                                sqrt((xPos[0] - xc1[0]) * (xPos[0] - xc1[0]) +
////                                     (xPos[1] - xc1[1]) * (xPos[1] - xc1[1])) > 0.6* rCore
////                                &&
////                                sqrt((xPos[0] - xc2[0]) * (xPos[0] - xc2[0]) +
////                                     (xPos[1] - xc2[1]) * (xPos[1] - xc2[1])) > 0.6* rCore
//                        );
//                    };
//                    auto fourthRuleHFiberP= [bDist](const TPZVec<REAL> &xPos) {
//                        return (xPos[0]>0.95*bDist ||
//                                xPos[1]>0.95*bDist);
//                    };
                    refineRulesP.Resize(2);
                    refineRulesP[0]=firstRuleHFiberP;
                    refineRulesP[1]=secondRuleHFiberP;
//                    refineRulesP[2]=thirdRuleHFiberP;
//                    refineRulesP[3]=fourthRuleHFiberP;
                }
                if(refineH==true){

                    auto firstRuleHFiberH = [xc1,xc2,rCore,bDist, Lambda](const TPZVec<REAL> &xPos) {
                        return (
                                sqrt(xPos[0]*xPos[0]+xPos[1]*xPos[1])<(Lambda + rCore)
                                &&
                                sqrt((xPos[0] - xc1[0]) * (xPos[0] - xc1[0]) +
                                     (xPos[1] - xc1[1]) * (xPos[1] - xc1[1])) > 0.67* rCore
                                &&
                                sqrt((xPos[0] - xc2[0]) * (xPos[0] - xc2[0]) +
                                     (xPos[1] - xc2[1]) * (xPos[1] - xc2[1])) > 0.67* rCore
                        );
                    };
//                    auto secondRuleHFiberH = [xc1,xc2,rCore,bDist](const TPZVec<REAL> &xPos) {
//                        return (
//                                sqrt(xPos[0]*xPos[0]+xPos[1]*xPos[1])<0.8*bDist
//                                &&
//                                sqrt((xPos[0] - xc1[0]) * (xPos[0] - xc1[0]) +
//                                     (xPos[1] - xc1[1]) * (xPos[1] - xc1[1])) > 0.8* rCore
//                                &&
//                                sqrt((xPos[0] - xc2[0]) * (xPos[0] - xc2[0]) +
//                                     (xPos[1] - xc2[1]) * (xPos[1] - xc2[1])) > 0.8* rCore
//                        );
//                    };
                    auto thirdRuleHFiberH= [bDist](const TPZVec<REAL> &xPos) {
                        return (xPos[0]>0.95*bDist ||
                                xPos[1]>0.95*bDist);
                    };
                    refineRulesH.Resize(2);
                    refineRulesH[0]=firstRuleHFiberH;
//                    refineRulesH[1]=secondRuleHFiberH;
//                    refineRulesH[2]=thirdRuleHFiberH;
                    refineRulesH[1]=thirdRuleHFiberH;
                }
                CreateHoleyFiberMesh(gmesh, simData.pzOpts.meshFile, matIdVec, simData.pzOpts.prefix,
                                     simData.pzOpts.exportGMesh, simData.pzOpts.scaleFactor,
                                     factor, boundTypeVec, pmlTypeVec, dPML, boundDist, symmetryX, symmetryY,
                                     refineH, refineRulesH);
                break;
            }
            case SPZModalAnalysisData::StepFiber:{
                const REAL realRCore = simData.physicalOpts.stepFiberOpts.realRCore;
                const REAL dPML = simData.physicalOpts.stepFiberOpts.dPML;
                const REAL boundDist = simData.physicalOpts.stepFiberOpts.boundDist;
                const REAL outerMaterialReffIndex = std::sqrt(std::real(simData.physicalOpts.erVec[1]));
                const int nLayersPML = simData.physicalOpts.stepFiberOpts.nLayersPml;
                CreateStepFiberMesh(gmesh, simData.pzOpts.meshFile, matIdVec, simData.pzOpts.prefix,
                                    simData.pzOpts.exportGMesh, simData.pzOpts.scaleFactor,
                                    factor, boundTypeVec, pmlTypeVec, realRCore, dPML,
                                    boundDist, outerMaterialReffIndex, nLayersPML);
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
    bool reallyRefineP = false;
    if(refineP == true) {
        if (simData.pzOpts.usingNeoPzMesh) {
            switch (simData.pzOpts.pzCase) {
                case SPZModalAnalysisData::HoleyFiber:{
                    reallyRefineP = true;
                    break;
                }
                default:
                    reallyRefineP = false;
                    break;
            }
        }
    }
    CreateCMesh(meshVec, gmesh, simData.pzOpts.pOrder, matIdVec, simData.physicalOpts.urVec, simData.physicalOpts.erVec,
                simData.physicalOpts.lambda, simData.physicalOpts.isCutOff, simData.pzOpts.prefix,
                simData.pzOpts.exportCMesh, simData.pzOpts.scaleFactor, simData.physicalOpts.alphaMax, boundTypeVec,
                pmlTypeVec, reallyRefineP, refineRulesP, simData.pzOpts.elType); // funcao para criar a malha computacional
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
    an.SetSolver(solver);

    if(simData.pzOpts.exportGMesh){
        std::cout<<"Printing mesh..."<<std::endl;
        TPZStack<std::string> scalnames, vecnames;
        vecnames.Push("Material");
        if(refineP){
            scalnames.Push("POrderH1");
            scalnames.Push("POrderHCurl");
        }
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
            const STATE currentKz = eigenValues[iSol];
            for(int j = 0; j < eigenVectors.Rows(); j++){
                currentEigenvector(j,0) = eigenVectors.GetVal(j,iSol);
            }
            for(int iMat = 0; iMat < matIdVec.size(); iMat++){
                TPZMatModalAnalysis *matPointer =
                        dynamic_cast<TPZMatModalAnalysis *>(cmesh->FindMaterial(matIdVec[iMat]));
                if(!matPointer) continue;
                matPointer->SetKz(currentKz);
                matPointer->SetPrintFieldRealPart(!simData.pzOpts.absVal);
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
        for(int icon = 0; icon < meshVec[k]->NConnects(); icon++){
            TPZConnect &con = meshVec[k]->ConnectVec()[icon];
            con.RemoveDepend();
        }
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

struct EdgeData {
    int coord1;
    int coord2;
    bool isNonLinearMapped;
    bool amIboundaryEdge;
    int quad1;
    int side1;
    int quad2;
    int side2;
};

//struct used for structural meshing of a general quadrilateral. the coordinates can then be used to create
//triangular or quadrilateral elements.
struct QuadrilateralData {
    TPZVec<int> pts;
    TPZVec<int> sides;
    TPZVec<bool> invSide;

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

    QuadrilateralData(const int &p1, const int &p2, const int &p3, const int &p4,
                      const TPZVec<REAL> &co1,const TPZVec<REAL> &co2,const TPZVec<REAL> &co3,const TPZVec<REAL> &co4) :
            coord1(co1),coord2(co2),coord3(co3),coord4(co4){
        auto map_quad_side_linear = [](const TPZVec<REAL> &coord1 ,const TPZVec<REAL> &coord2,const REAL & s) {
            TPZVec<REAL> point(2,0.);
            point[0] = (coord1[0] - s*coord1[0] + coord2[0] + s*coord2[0])/2.;
            point[1] = (coord1[1] - s*coord1[1] + coord2[1] + s*coord2[1])/2.;
            return point;
        };
        pts.Resize(4);
        pts[0]=p1;
        pts[1]=p2;
        pts[2]=p3;
        pts[3]=p4;

        isSide1nonLinear = false;
        mapSide1 = std::bind(map_quad_side_linear,coord1,coord2,std::placeholders::_1);
        isSide2nonLinear = false;
        mapSide2 = std::bind(map_quad_side_linear,coord2,coord3,std::placeholders::_1);
        isSide3nonLinear = false;
        mapSide3 = std::bind(map_quad_side_linear,coord3,coord4,std::placeholders::_1);
        isSide4nonLinear = false;
        mapSide4 = std::bind(map_quad_side_linear,coord4,coord1,std::placeholders::_1);

    }

    QuadrilateralData(const TPZVec<REAL> &p1, const TPZVec<REAL> &p2, const TPZVec<REAL> &p3, const TPZVec<REAL> &p4) :
            coord1(p1),coord2(p2),coord3(p3),coord4(p4){
        auto map_quad_side_linear = [](const TPZVec<REAL> &coord1 ,const TPZVec<REAL> &coord2,const REAL & s) {
            TPZVec<REAL> point(2,0.);
            point[0] = (coord1[0] - s*coord1[0] + coord2[0] + s*coord2[0])/2.;
            point[1] = (coord1[1] - s*coord1[1] + coord2[1] + s*coord2[1])/2.;
            return point;
        };

        pts.Resize(0);
        sides.Resize(0);
        invSide.Resize(0);

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

    void SetMapsLinearity(const bool &m1, const bool &m2, const bool &m3, const bool &m4, const int quad, TPZVec<EdgeData> &allEdges){
        isSide1nonLinear = m1;
        isSide2nonLinear = m2;
        isSide3nonLinear = m3;
        isSide4nonLinear = m4;
        sides.Resize(4);
        invSide.Resize(4);
        for(int iSide = 0; iSide < 4; iSide++){
            const int &pt1 = pts[iSide];
            const int &pt2 = pts[(iSide+1)%4];
            int foundEdge = -1;
            for(int iEdge = 0; iEdge < allEdges.size(); iEdge++){
                if(allEdges[iEdge].coord1 == pt1 && allEdges[iEdge].coord2 == pt2) {
                    std::cout<<"found existing edge with same orientation"<<std::endl;
                    DebugStop();
                }
                if(allEdges[iEdge].coord1 == pt2 && allEdges[iEdge].coord2 == pt1) {
                    foundEdge = iEdge;
                    invSide = true;
                    allEdges[foundEdge].quad2 = quad;
                    allEdges[foundEdge].side2 = iSide + 1;
                    allEdges[foundEdge].amIboundaryEdge = false;
                }
            }
            if(foundEdge == -1){
                foundEdge = allEdges.size();
                invSide[iSide] = false;
                allEdges.Resize(foundEdge + 1);
                allEdges[foundEdge].coord1 = pt1;
                allEdges[foundEdge].coord2 = pt2;
                allEdges[foundEdge].quad1 = quad;
                allEdges[foundEdge].side1 = iSide + 1;
                allEdges[foundEdge].quad2 = quad;
                allEdges[foundEdge].side2 = iSide + 1;
                allEdges[foundEdge].amIboundaryEdge = true;
                switch(iSide){
                    case 0:
                        allEdges[foundEdge].isNonLinearMapped = isSide1nonLinear;
                        break;
                    case 1:
                        allEdges[foundEdge].isNonLinearMapped = isSide2nonLinear;
                        break;
                    case 2:
                        allEdges[foundEdge].isNonLinearMapped = isSide3nonLinear;
                        break;
                    case 3:
                        allEdges[foundEdge].isNonLinearMapped = isSide4nonLinear;
                        break;
                    default:
                        DebugStop();
                }
            }
            sides[iSide] = foundEdge;
        }
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
                                     const SPZModalAnalysisData::boundtype &symmetryType, bool refine) {
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
    if(refine){
        const REAL margin = 0.2 * wDomain * scale;
        TPZVec<REAL> qsiPos(2,0.25);
        TPZVec<REAL> xPos;
        TPZVec<TPZGeoEl *>sons;
        int nel = gmesh->NElements();
        for(int iel =0 ; iel< nel; iel++){
            TPZGeoEl *geo = gmesh->Element(iel);
            geo->X(qsiPos,xPos);
            if(std::abs(xPos[0]-wDomain*scale/2) < margin){
                if(!geo->Father()){
                    geo->Divide(sons);
                    nel = gmesh->NElements();
                }
            }
        }
    }
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
                    const REAL &boundDist, const REAL &outerReffIndex, const int &nLayersPml, bool refine) {
    const int nDivTCore = factor * 2, nDivRCore = factor * 5, nDivTCladding = factor * 4, nDivPml = factor * nLayersPml + 1;

    if(std::min<int>({nDivTCore,nDivRCore,nDivTCladding,nDivPml}) < 2 ) {
        std::cout<<"Mesh has not sufficient divisions."<<std::endl;
        std::cout<<"Minimum is 2."<<std::endl;
        DebugStop();
    }

    const int nQuads = 17;
    const int nPoints = 24;

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
    const REAL bound = rCore + boundDist * outerReffIndex * 2 * M_PI;
    TPZManVector<REAL,2> xc(2, 0.);
    xc[0] = 0.;
    xc[1] = 0.;



    ///all points
    TPZManVector<TPZVec<REAL>,33> pointsVec(nPoints,TPZVec<REAL>(2,0.));
    TPZManVector<EdgeData,33> edgesVec(0,EdgeData());
    //hole 1
    pointsVec[0][0]= xc[0] + M_SQRT1_2 * std::cos(M_PI_4)*rCore; pointsVec[0][1]= xc[1] + M_SQRT1_2 * std::sin(M_PI_4)*rCore;
    pointsVec[1][0]= xc[0] - M_SQRT1_2 * std::cos(M_PI_4)*rCore; pointsVec[1][1]= xc[1] + M_SQRT1_2 * std::sin(M_PI_4)*rCore;
    pointsVec[2][0]= xc[0] - M_SQRT1_2 * std::cos(M_PI_4)*rCore; pointsVec[2][1]= xc[1] - M_SQRT1_2 * std::sin(M_PI_4)*rCore;
    pointsVec[3][0]= xc[0] + M_SQRT1_2 * std::cos(M_PI_4)*rCore; pointsVec[3][1]= xc[1] - M_SQRT1_2 * std::sin(M_PI_4)*rCore;

    pointsVec[4][0]= xc[0] + std::cos(M_PI_4)*rCore; pointsVec[4][1]= xc[1] + std::sin(M_PI_4)*rCore;
    pointsVec[5][0]= xc[0] - std::cos(M_PI_4)*rCore; pointsVec[5][1]= xc[1] + std::sin(M_PI_4)*rCore;
    pointsVec[6][0]= xc[0] - std::cos(M_PI_4)*rCore; pointsVec[6][1]= xc[1] - std::sin(M_PI_4)*rCore;
    pointsVec[7][0]= xc[0] + std::cos(M_PI_4)*rCore; pointsVec[7][1]= xc[1] - std::sin(M_PI_4)*rCore;


    //cladding
    pointsVec[8][0] =  1 * bound; pointsVec[8][1] = -1 * bound;
    pointsVec[9][0] =  1 * bound; pointsVec[9][1] =  1 * bound;
    pointsVec[10][0] = -1 * bound; pointsVec[10][1] =  1 * bound;
    pointsVec[11][0] = -1 * bound; pointsVec[11][1] = -1 * bound;
    //PML
    pointsVec[12][0] =  1 * (bound+lengthPML); pointsVec[12][1] = -1 * bound;
    pointsVec[13][0] =  1 * (bound+lengthPML); pointsVec[13][1] =  1 * bound;

    pointsVec[14][0] =  1 * bound; pointsVec[14][1] =  1 * (bound+lengthPML);
    pointsVec[15][0] = -1 * bound; pointsVec[15][1] =  1 * (bound+lengthPML);

    pointsVec[16][0] = -1 * (bound+lengthPML); pointsVec[16][1] =  1 * bound;
    pointsVec[17][0] = -1 * (bound+lengthPML); pointsVec[17][1] = -1 * bound;

    pointsVec[18][0] = -1 * bound; pointsVec[18][1] = -1 * (bound+lengthPML);
    pointsVec[19][0] = 1 * bound; pointsVec[19][1] = -1 * (bound+lengthPML);

    pointsVec[20][0] =  1 * (bound+lengthPML); pointsVec[20][1] = -1 * (bound+lengthPML);
    pointsVec[21][0] =  1 * (bound+lengthPML); pointsVec[21][1] =  1 * (bound+lengthPML);
    pointsVec[22][0] = -1 * (bound+lengthPML); pointsVec[22][1] =  1 * (bound+lengthPML);
    pointsVec[23][0] = -1 * (bound+lengthPML); pointsVec[23][1] = -1 * (bound+lengthPML);

    TPZFMatrix<int> quadPointsVec(nQuads,4);
    quadPointsVec(0,0) = 0; quadPointsVec(0,1) = 1; quadPointsVec(0,2) = 2; quadPointsVec(0,3) = 3;
    quadPointsVec(1,0) = 4; quadPointsVec(1,1) = 0; quadPointsVec(1,2) = 3; quadPointsVec(1,3) = 7;
    quadPointsVec(2,0) = 4; quadPointsVec(2,1) = 5; quadPointsVec(2,2) = 1; quadPointsVec(2,3) = 0;
    quadPointsVec(3,0) = 1; quadPointsVec(3,1) = 5; quadPointsVec(3,2) = 6; quadPointsVec(3,3) = 2;
    quadPointsVec(4,0) = 3; quadPointsVec(4,1) = 2; quadPointsVec(4,2) = 6; quadPointsVec(4,3) = 7;

    quadPointsVec(5,0) = 9; quadPointsVec(5,1) = 4; quadPointsVec(5,2) = 7; quadPointsVec(5,3) = 8;
    quadPointsVec(6,0) = 9; quadPointsVec(6,1) = 10; quadPointsVec(6,2) = 5; quadPointsVec(6,3) = 4;
    quadPointsVec(7,0) = 5; quadPointsVec(7,1) = 10; quadPointsVec(7,2) = 11; quadPointsVec(7,3) = 6;
    quadPointsVec(8,0) = 7; quadPointsVec(8,1) = 6; quadPointsVec(8,2) = 11; quadPointsVec(8,3) = 8;

    quadPointsVec(9,0) = 13; quadPointsVec(9,1) = 9; quadPointsVec(9,2) = 8; quadPointsVec(9,3) = 12;
    quadPointsVec(10,0) = 14; quadPointsVec(10,1) = 15; quadPointsVec(10,2) = 10; quadPointsVec(10,3) = 9;
    quadPointsVec(11,0) = 10; quadPointsVec(11,1) = 16; quadPointsVec(11,2) = 17; quadPointsVec(11,3) = 11;
    quadPointsVec(12,0) = 8; quadPointsVec(12,1) = 11; quadPointsVec(12,2) = 18; quadPointsVec(12,3) = 19;

    quadPointsVec(13,0) = 21; quadPointsVec(13,1) = 14; quadPointsVec(13,2) =  9; quadPointsVec(13,3) = 13;
    quadPointsVec(14,0) = 15; quadPointsVec(14,1) = 22; quadPointsVec(14,2) = 16; quadPointsVec(14,3) = 10;
    quadPointsVec(15,0) = 11; quadPointsVec(15,1) = 17; quadPointsVec(15,2) = 23; quadPointsVec(15,3) = 18;
    quadPointsVec(16,0) = 12; quadPointsVec(16,1) =  8; quadPointsVec(16,2) = 19; quadPointsVec(16,3) = 20;

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
    //first hole
    nDivQsi[0]=nDivRCore;nDivEta[0]=nDivRCore;
    nDivQsi[1]=nDivTCore; nDivEta[1]=nDivRCore;
    nDivQsi[2]=nDivRCore; nDivEta[2]=nDivTCore;
    nDivQsi[3]=nDivTCore; nDivEta[3]=nDivRCore;
    nDivQsi[4]=nDivRCore; nDivEta[4]=nDivTCore;
    //cladding
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

    side4NonLinearVec[1] = true;
    side1NonLinearVec[2] = true;
    side2NonLinearVec[3] = true;
    side3NonLinearVec[4] = true;

    side2NonLinearVec[5] = true;
    side3NonLinearVec[6] = true;
    side4NonLinearVec[7] = true;
    side1NonLinearVec[8] = true;


    TPZManVector<TPZVec<REAL>,14> thetaVec(8,TPZVec<REAL>(2,0.));

    thetaVec[0][0] =-1 * M_PI/4;thetaVec[0][1] = 1 * M_PI/4;//quad1
    thetaVec[1][0] = 1 * M_PI/4;thetaVec[1][1] = 3 * M_PI/4;//quad2
    thetaVec[2][0] = 3 * M_PI/4;thetaVec[2][1] = 5 * M_PI/4;//quad3
    thetaVec[3][0] = 5 * M_PI/4;thetaVec[3][1] = 7 * M_PI/4;//quad4

    thetaVec[4][0] = 1 * M_PI/4;thetaVec[4][1] =-1 * M_PI/4;//quad5
    thetaVec[5][0] = 3 * M_PI/4;thetaVec[5][1] = 1 * M_PI/4;//quad6
    thetaVec[6][0] = 5 * M_PI/4;thetaVec[6][1] = 3 * M_PI/4;//quad7
    thetaVec[7][0] = 7 * M_PI/4;thetaVec[7][1] = 5 * M_PI/4;//quad8


    TPZVec<TPZVec<REAL>*> xcRef(8,&xc);

    TPZVec<int> matIdBoundVec(4,-1);
    matIdBoundVec[0] = matIdBC;//low
    matIdBoundVec[1] = matIdBC;//down
    matIdBoundVec[2] = matIdBC;//up
    matIdBoundVec[3] = matIdBC;//left

    TPZVec<REAL> boundDistVec(4,-1);
    boundDistVec[0] = 0.;//low
    boundDistVec[1] = bound+lengthPML;//down
    boundDistVec[2] = bound+lengthPML;//up
    boundDistVec[3] = 0.;//left

    TPZVec<REAL> rVec(14,rCore);

    gmesh = CreateStructuredMesh(pointsVec, edgesVec, quadPointsVec, matIdsQuads, nDivQsi,
                                 nDivEta, side1NonLinearVec, side2NonLinearVec, side3NonLinearVec, side4NonLinearVec,
                                 thetaVec, xcRef, matIdBoundVec, boundDistVec, rVec);
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

    if(refine){
        TPZRefPatternDataBase db;
        db.InitializeUniformRefPattern(EQuadrilateral);
        const REAL margin = 0.4 * rCore;
        TPZVec<REAL> qsiPos(2,0);
        TPZVec<REAL> xPos;
        TPZVec<TPZGeoEl *>sons;
        int nel = gmesh->NElements();
        for(int iel =0 ; iel< nel; iel++){
            TPZGeoElRefPattern<pzgeom::TPZGeoQuad> *geo = dynamic_cast<TPZGeoElRefPattern<pzgeom::TPZGeoQuad> *> (gmesh->Element(iel));
            if(!geo) continue;
            geo->X(qsiPos,xPos);
            const REAL dist = sqrt((xPos[0]-xc[0])*(xPos[0]-xc[0])+(xPos[1]-xc[1])*(xPos[1]-xc[1]));
            if(std::abs(rCore - dist) < margin){
                if(geo->Type() == EQuadrilateral && geo->HasSubElement() == 0){
                    if(!geo->Father()){
                        auto oldref = geo->GetRefPattern();
                        auto uniformref = db.GetUniformRefPattern(EQuadrilateral);
                        geo->SetRefPattern(uniformref);
                        geo->Divide(sons);
                        //geo->SetRefPattern(oldref);
                        for(int i = 0; i < sons.size(); i++){
                            sons[i]->SetRefPattern(oldref);
                        }
                        nel = gmesh->NElements();
                    }
                }
            }
        }
    }
    gmesh->BuildConnectivity();

    for (long iel = 0; iel < nel; iel++) {
        TPZGeoEl *gelQuad = gmesh->ElementVec()[iel];
        if(gelQuad->Type() == EQuadrilateral && gelQuad->HasSubElement() == 0){
            gelQuad->Divide(sons); nel = gmesh->NElements(); continue;
        }
    }

    gmesh->BuildConnectivity();

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

TPZGeoMesh *
CreateStructuredMesh(const TPZVec<TPZVec<REAL>> &pointsVec, TPZVec<EdgeData> &edgesVec, const TPZFMatrix<int> &quadPointsVec,
                     const TPZVec<int> &matIdsQuads, const TPZVec<int> &nDivQsi, const TPZVec<int> &nDivEta,
                     const TPZVec<bool> &side1NonLinearVec, const TPZVec<bool> &side2NonLinearVec,
                     const TPZVec<bool> &side3NonLinearVec, const TPZVec<bool> &side4NonLinearVec,
                     const TPZVec<TPZVec<REAL>> &thetaVec, const TPZVec<TPZVec<REAL> *> &xcRef,
                     const TPZVec<int> &matIdBoundVec, const TPZVec<REAL> &boundDistVec, const TPZVec<REAL> &rVec,
                     const bool nonLinearMapping) {

        auto map_quad_side_arc = [](const TPZVec<REAL> &theta ,const TPZVec<REAL> &xc, const REAL &r, const REAL & s) {
            TPZVec<REAL> point(2,0.);
            point[0] = xc[0] + r*cos((theta[1] + s*theta[1] + theta[0] - s*theta[0])/2.);
            point[1] = xc[1] + r*sin((theta[1] + s*theta[1] + theta[0] - s*theta[0])/2.);
            return point;
        };
        TPZGeoMesh *gmesh;
        const int nQuads = side1NonLinearVec.size();

        TPZVec<QuadrilateralData *> quadVec(nQuads, nullptr);
        int iNonLinear = 0;
        for (int iQuad = 0; iQuad < nQuads; iQuad++) {
            quadVec[iQuad] = new QuadrilateralData(quadPointsVec.GetVal(iQuad, 0), quadPointsVec.GetVal(iQuad, 1),
                                                   quadPointsVec.GetVal(iQuad, 2), quadPointsVec.GetVal(iQuad, 3),
                                                   pointsVec[quadPointsVec.GetVal(iQuad, 0)],
                                                   pointsVec[quadPointsVec.GetVal(iQuad, 1)],
                                                   pointsVec[quadPointsVec.GetVal(iQuad, 2)],
                                                   pointsVec[quadPointsVec.GetVal(iQuad, 3)]);
            if (side1NonLinearVec[iQuad]) {
                TPZVec<REAL> &xc = *(xcRef[iNonLinear]);
                quadVec[iQuad]->mapSide1 = std::__1::bind(map_quad_side_arc, thetaVec[iNonLinear], xc, rVec[iNonLinear],
                                                          std::__1::placeholders::_1);
                iNonLinear++;
            } else if (side2NonLinearVec[iQuad]) {
                TPZVec<REAL> &xc = *(xcRef[iNonLinear]);
                quadVec[iQuad]->mapSide2 = std::__1::bind(map_quad_side_arc, thetaVec[iNonLinear], xc, rVec[iNonLinear],
                                                          std::__1::placeholders::_1);
                iNonLinear++;
            } else if (side3NonLinearVec[iQuad]) {
                TPZVec<REAL> &xc = *(xcRef[iNonLinear]);
                quadVec[iQuad]->mapSide3 = std::__1::bind(map_quad_side_arc, thetaVec[iNonLinear], xc, rVec[iNonLinear],
                                                          std::__1::placeholders::_1);
                iNonLinear++;
            } else if (side4NonLinearVec[iQuad]) {
                TPZVec<REAL> &xc = *(xcRef[iNonLinear]);
                quadVec[iQuad]->mapSide4 = std::__1::bind(map_quad_side_arc, thetaVec[iNonLinear], xc, rVec[iNonLinear],
                                                          std::__1::placeholders::_1);
                iNonLinear++;
            }
            quadVec[iQuad]->SetMapsLinearity(side1NonLinearVec[iQuad], side2NonLinearVec[iQuad],
                                             side3NonLinearVec[iQuad], side4NonLinearVec[iQuad],
                                             iQuad, edgesVec);
            quadVec[iQuad]->CreateQuadrilateral(nDivQsi[iQuad], nDivEta[iQuad]);
        }

        ////////////////////////////////////////CREATE NODES///////////////////////////////////////
        ////////////////////////////////////////FOR INTERIOR///////////////////////////////////////
        ///////////////////////////////////////////POINTS//////////////////////////////////////////
        long nodeId = 0;
        gmesh = new TPZGeoMesh;
        TPZVec<TPZVec<long>> elNodeFaceVecs(nQuads, TPZVec<long>(0, -1));

        for (int el = 0; el < nQuads; el++) {
            const long nNodesOld = gmesh->NodeVec().NElements();
            const long nNodesEl = (nDivQsi[el] - 2) * (nDivEta[el] - 2);
            gmesh->NodeVec().Resize(nNodesOld + nNodesEl);
            elNodeFaceVecs[el].Resize(nDivEta[el] * nDivQsi[el], -1);
            //interior nodes only
            for (int i = 1; i < nDivEta[el] - 1; i++) {
                for (int j = 1; j < nDivQsi[el] - 1; j++) {
                    TPZManVector<REAL, 3> node(3, 0.);
                    node[0] = quadVec[el]->facePts(i * nDivQsi[el] + j, 0);
                    node[1] = quadVec[el]->facePts(i * nDivQsi[el] + j, 1);
                    gmesh->NodeVec()[nodeId].Initialize(node, *gmesh);
                    elNodeFaceVecs[el][i * nDivQsi[el] + j] = nodeId;
                    nodeId++;
                }
            }
        }

//    std::cout<<"NUMBER OF EDGES "<<edgesVec.size()<<std::endl;
//    for(int iEdge = 0; iEdge < edgesVec.size(); iEdge++){
//        std::cout<<"Edge "<<iEdge<<std::setfill(' ');;
//        std::cout<<"\tquad1:"<<std::setw(2)<<edgesVec[iEdge].quad1<<"\tquad2:"<<std::setw(2)<<edgesVec[iEdge].quad2;
//        std::cout<<"\tNon-linear mapped:"<<edgesVec[iEdge].isNonLinearMapped;
//        std::cout<<"\tBoundary:"<<edgesVec[iEdge].amIboundaryEdge<<std::endl;
//    }
        const int nEdges = edgesVec.size();
        ////////////////////////////////////////CREATE NODES///////////////////////////////////////
        //////////////////////////////////////////FOR EDGE/////////////////////////////////////////
        ///////////////////////////////////////////POINTS//////////////////////////////////////////
        auto sidePos = [](const int side, const int &i, const int &nQsi, const int &nEta, const bool &reverse) {
            const int iP = reverse ? (side % 2 ? nQsi - 1 - i : nEta - 1 - i) : i;
            switch (side) {
                case 1 :
                    return iP;
                case 2 :
                    return (iP * nQsi + (nQsi - 1));
                case 3 :
                    return (nQsi - 1 - iP + nQsi * (nEta - 1));
                case 4 :
                    return (nEta - 1 - iP) * nQsi;
                default:
                    DebugStop();
                    return -1;
            }
        };

        auto createEdgeNodes = [sidePos, nDivEta, nDivQsi, quadVec, gmesh]
                (TPZVec<TPZVec<long>> &elNodeFaceVecs, const long &el1, const long &el2,
                 const int &side1, const int &side2, long &nodeId,
                 TPZVec<long> &side1pts, TPZVec<long> &side2pts, const bool &revEl2) {
            const int nPts = side1 % 2 ? nDivQsi[el1] : nDivEta[el1];
            long nNodesOld = gmesh->NodeVec().NElements();
            long nNodesEl = nPts;
            for (int i = 0; i < nPts; i++) {
                const int posEl1 = sidePos(side1, i, nDivQsi[el1], nDivEta[el1], false);
                const int posEl2 = sidePos(side2, i, nDivQsi[el2], nDivEta[el2], revEl2);
                if (elNodeFaceVecs[el1][posEl1] != -1) {
                    elNodeFaceVecs[el2][posEl2] = elNodeFaceVecs[el1][posEl1];
                    continue;
                }
                if (elNodeFaceVecs[el2][posEl2] != -1) {
                    elNodeFaceVecs[el1][posEl1] = elNodeFaceVecs[el2][posEl2];
                    continue;
                }
                TPZManVector<REAL, 3> node(3, 0.);
                node[0] = quadVec[el1]->facePts(posEl1, 0);
                node[1] = quadVec[el1]->facePts(posEl1, 1);
                long nNodesOld = gmesh->NodeVec().NElements();
                gmesh->NodeVec().Resize(nNodesOld + 1);
                gmesh->NodeVec()[nodeId].Initialize(node, *gmesh);
                elNodeFaceVecs[el1][posEl1] = nodeId;
                elNodeFaceVecs[el2][posEl2] = nodeId;
                nodeId++;
            }
            TPZFMatrix<REAL> *nonLinearPts = nullptr;
            if (quadVec[el1]->isSide1nonLinear && side1 == 1) { nonLinearPts = &quadVec[el1]->side1IntPts; }
            else if (quadVec[el1]->isSide2nonLinear && side1 == 2) { nonLinearPts = &quadVec[el1]->side2IntPts; }
            else if (quadVec[el1]->isSide3nonLinear && side1 == 3) { nonLinearPts = &quadVec[el1]->side3IntPts; }
            else if (quadVec[el1]->isSide4nonLinear && side1 == 4) { nonLinearPts = &quadVec[el1]->side4IntPts; }
            long nNonLinPts = 0;
            if (nonLinearPts) {
                nNonLinPts = nonLinearPts->Rows();
                side1pts.Resize(nNonLinPts);
                side2pts.Resize(nNonLinPts);
                nNodesOld = gmesh->NodeVec().NElements();
                nNodesEl = nNonLinPts;
                gmesh->NodeVec().Resize(nNodesOld + nNodesEl);
            }
            for (int i = 0; i < nNonLinPts; i++) {
                const int posEl1 = i;
                const int posEl2 = revEl2 ? nNonLinPts - 1 - i : i;
                TPZManVector<REAL, 3> node(3, 0.);
                node[0] += (*nonLinearPts)(posEl1, 0);
                node[1] += (*nonLinearPts)(posEl1, 1);
                gmesh->NodeVec()[nodeId].Initialize(node, *gmesh);
                side1pts[posEl1] = nodeId;
                side2pts[posEl2] = nodeId;
                nodeId++;
            }
        };

        TPZVec<TPZVec<long>> elNodeSide1Vecs(nQuads, TPZVec<long>(0, -1));
        TPZVec<TPZVec<long>> elNodeSide2Vecs(nQuads, TPZVec<long>(0, -1));
        TPZVec<TPZVec<long>> elNodeSide3Vecs(nQuads, TPZVec<long>(0, -1));
        TPZVec<TPZVec<long>> elNodeSide4Vecs(nQuads, TPZVec<long>(0, -1));
        for (int i = 0; i < nEdges; i++) {
            TPZVec<long> *side1pts = nullptr;
            TPZVec<long> *side2pts = nullptr;
            int quad1 = edgesVec[i].quad1, quad2 = edgesVec[i].quad2,
                    side1 = edgesVec[i].side1, side2 = edgesVec[i].side2;
            bool revEl = !(edgesVec[i].amIboundaryEdge);

            switch (side1) {
                case 1:
                    side1pts = &elNodeSide1Vecs[quad1];
                    break;
                case 2:
                    side1pts = &elNodeSide2Vecs[quad1];
                    break;
                case 3:
                    side1pts = &elNodeSide3Vecs[quad1];
                    break;
                case 4:
                    side1pts = &elNodeSide4Vecs[quad1];
                    break;
            }

            switch (side2) {
                case 1:
                    side2pts = &elNodeSide1Vecs[quad2];
                    break;
                case 2:
                    side2pts = &elNodeSide2Vecs[quad2];
                    break;
                case 3:
                    side2pts = &elNodeSide3Vecs[quad2];
                    break;
                case 4:
                    side2pts = &elNodeSide4Vecs[quad2];
                    break;
            }
            createEdgeNodes(elNodeFaceVecs, quad1, quad2, side1, side2, nodeId,
                            *side1pts, *side2pts, revEl);
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
        std::__1::istringstream str(buf);
        refp = new TPZRefPattern(str);
        refp->GenerateSideRefPatterns();
        gRefDBase.InsertRefPattern(refp);
        gRefDBase.InitializeUniformRefPattern(EOned);
        if (!refp) {
            DebugStop();
        }
        for (int quad = 0; quad < nQuads; quad++) {
            for (int i = 0; i < nDivEta[quad] - 1; i++) {
                for (int j = 0; j < nDivQsi[quad] - 1; j++) {

                    const int node0 = i * nDivQsi[quad] + j;//Lower left vertex
                    const int node1 = i * nDivQsi[quad] + j + 1;//Lower right vertex
                    const int node2 = (i + 1) * nDivQsi[quad] + j + 1;//Upper right vertex
                    const int node3 = (i + 1) * nDivQsi[quad] + j;//Upper left vertex

                    TPZManVector<long, 4> nodesIdVec(4, -1);

                    const int matId = matIdsQuads[quad];
                    nodesIdVec[0] = elNodeFaceVecs[quad][node0];
                    nodesIdVec[1] = elNodeFaceVecs[quad][node1];
                    nodesIdVec[2] = elNodeFaceVecs[quad][node2];
                    nodesIdVec[3] = elNodeFaceVecs[quad][node3];
                    if (nodesIdVec[0] == -1 || nodesIdVec[1] == -1 || nodesIdVec[2] == -1 || nodesIdVec[3] == -1) {
                        DebugStop();
                    }
                    if ((i == 0 && quadVec[quad]->isSide1nonLinear) ||
                        (j == nDivQsi[quad] - 2 && quadVec[quad]->isSide2nonLinear) ||
                        (i == nDivEta[quad] - 2 && quadVec[quad]->isSide3nonLinear) ||
                        (j == 0 && quadVec[quad]->isSide4nonLinear)) {
                        TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> > *quadEl =
                                new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> >(nodesIdVec, matId,
                                                                                                 *gmesh);
                        quadEl->SetRefPattern(refp);
                    } else {
                        TPZGeoElRefPattern<pzgeom::TPZGeoQuad> *quadEl =
                                new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodesIdVec, matId, *gmesh);
                        quadEl->SetRefPattern(refp);
                    }
                }
            }
        }

        for (int edge = 0; edge < nEdges; edge++) {
            auto refpArc = gRefDBase.GetUniformRefPattern(EOned);
            int quad = edgesVec[edge].quad1, side = edgesVec[edge].side1;
            if (edgesVec[edge].amIboundaryEdge) {//boundary edge
                const int nArcs = side % 2 ? nDivQsi[quad] - 1 : nDivEta[quad] - 1;
                TPZVec<REAL> &pt1 = pointsVec[edgesVec[edge].coord1];
                TPZVec<REAL> &pt2 = pointsVec[edgesVec[edge].coord2];
                REAL tol = 1e-08;
                int matId = -666;
                if (std::__1::abs(pt1[1] - pt2[1]) < tol) {//horizontal edge
                    if (std::__1::abs(pt1[1] - boundDistVec[2]) < tol) {
                        matId = matIdBoundVec[2];
                    } else {
                        matId = matIdBoundVec[0];
                    }
                } else if (std::__1::abs(pt1[0] - pt2[0]) < tol) {//vertical edge
                    if (std::__1::abs(pt1[0] - boundDistVec[1]) < tol) {
                        matId = matIdBoundVec[1];
                    } else {
                        matId = matIdBoundVec[3];
                    }
                } else {
                    DebugStop();
                }
                for (int i = 0; i < nArcs; i++) {
                    TPZManVector<long, 3> nodesIdVec(2, 0.);
                    const int vertex1 = sidePos(side, i, nDivQsi[quad], nDivEta[quad], false);
                    const int vertex2 = sidePos(side, i + 1, nDivQsi[quad], nDivEta[quad], false);


                    nodesIdVec[0] = elNodeFaceVecs[quad][vertex1];
                    nodesIdVec[1] = elNodeFaceVecs[quad][vertex2];
                    if (nodesIdVec[0] == -1 || nodesIdVec[1] == -1) {
                        DebugStop();
                    }

                    TPZGeoElRefPattern<pzgeom::TPZGeoLinear> *arc =
                            new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesIdVec, matId, *gmesh);
                    arc->SetRefPattern(refpArc);
                }
                continue;
            }

            TPZVec<long> *intPointsCoord = nullptr;
            switch (side) {
                case 1:
                    if (quadVec[quad]->isSide1nonLinear == false) continue;
                    intPointsCoord = &(elNodeSide1Vecs[quad]);
                    break;
                case 2:
                    if (quadVec[quad]->isSide2nonLinear == false) continue;
                    intPointsCoord = &(elNodeSide2Vecs[quad]);
                    break;
                case 3:
                    if (quadVec[quad]->isSide3nonLinear == false) continue;
                    intPointsCoord = &(elNodeSide3Vecs[quad]);
                    break;
                case 4:
                    if (quadVec[quad]->isSide4nonLinear == false) continue;
                    intPointsCoord = &(elNodeSide4Vecs[quad]);
                    break;
            }

            const int nArcs = side % 2 ? nDivQsi[quad] - 1 : nDivEta[quad] - 1;
            //auto sidePos = [](const int side, const int &i, const int &nQsi, const int &nEta, const bool &reverse){
            if(nonLinearMapping == false) continue;
            for (int i = 0; i < nArcs; i++) {
                TPZManVector<long, 3> nodesIdVec(3, 0.);
                const int vertex1 = sidePos(side, i, nDivQsi[quad], nDivEta[quad], false);
                const int vertex2 = sidePos(side, i + 1, nDivQsi[quad], nDivEta[quad], false);


                nodesIdVec[0] = elNodeFaceVecs[quad][vertex1];
                nodesIdVec[1] = elNodeFaceVecs[quad][vertex2];
                nodesIdVec[2] = (*intPointsCoord)[i];//mid-point
                if (nodesIdVec[0] == -1 || nodesIdVec[1] == -1 || nodesIdVec[2] == -1) {
                    DebugStop();
                }
                TPZGeoElRefPattern<pzgeom::TPZArc3D> *arc =
                        new TPZGeoElRefPattern<pzgeom::TPZArc3D>(nodesIdVec, -24, *gmesh);//random material id
                arc->SetRefPattern(refpArc);
            }
        }

    return gmesh;
    }

void
CreateHoleyFiberMesh(TPZGeoMesh *&gmesh, const std::string mshFileName, TPZVec<int> &matIdVec,
                     const std::string &prefix, const bool &print, const REAL &scale, const int &factor,
                     TPZVec<SPZModalAnalysisData::boundtype> &boundTypeVec,
                     TPZVec<SPZModalAnalysisData::pmltype> &pmlTypeVec, const REAL &dPML, const REAL &boundDist,
                     SPZModalAnalysisData::boundtype symmetryX, SPZModalAnalysisData::boundtype symmetryY, bool refine,
                     TPZVec<std::function<bool (const TPZVec<REAL> &)>> refineRules) {
    //SPZModalAnalysisData::boundtype symmetry;
    const int nDivTHole = factor * 2 + 1, nDivRHole = factor * 4 + 1,
            nDivCladding1 = factor * 6 + 1, nDivCladding2 = factor * 4+1,nDivCladding3 = factor * 2+1,
                    nDivPml = factor * 3 + 1;

    const int nQuads = 32;
    const int nPoints = 44;

    if(std::min<int>({nDivTHole,nDivRHole,nDivCladding1,nDivPml}) < 2 ) {
        std::cout<<"Mesh has not sufficient divisions."<<std::endl;
        std::cout<<"Minimum is 2."<<std::endl;
        DebugStop();
    }

    const int matIdHole = 1, matIdCladding = 2, matIdBC= 8,matIdSymmetryX = 7, matIdSymmetryY = 6;
    const int matIdPMLxp = 3,
            matIdPMLyp = 4,
            matIdPMLxpyp = 5;

    const REAL rCore = scale * 2.5e-6;
    const REAL lengthPML = scale * dPML;
    const REAL bound = scale * boundDist;
    const REAL Lambda = scale * 6.75e-6;
    //holes centers
    TPZManVector<REAL,2> xc1(2, 0.);
    TPZManVector<REAL,2> xc2(2, 0.);
    xc1[0] = Lambda * 0.5;
    xc1[1] = Lambda * std::sqrt(3) * 0.5;
    xc2[0] = Lambda;
    xc2[1] = 0;



    ///all points
    TPZManVector<TPZVec<REAL>,33> pointsVec(nPoints,TPZVec<REAL>(2,0.));
    TPZManVector<EdgeData,33> edgesVec(0,EdgeData());
    //hole 1
    pointsVec[0][0]= xc1[0] + 0.5 * std::cos(M_PI_4)*rCore; pointsVec[0][1]= xc1[1] + 0.5 * std::sin(M_PI_4)*rCore;
    pointsVec[1][0]= xc1[0] - 0.5 * std::cos(M_PI_4)*rCore; pointsVec[1][1]= xc1[1] + 0.5 * std::sin(M_PI_4)*rCore;
    pointsVec[2][0]= xc1[0] - 0.5 * std::cos(M_PI_4)*rCore; pointsVec[2][1]= xc1[1] - 0.5 * std::sin(M_PI_4)*rCore;
    pointsVec[3][0]= xc1[0] + 0.5 * std::cos(M_PI_4)*rCore; pointsVec[3][1]= xc1[1] - 0.5 * std::sin(M_PI_4)*rCore;

    pointsVec[4][0]= xc1[0] + std::cos(M_PI_4)*rCore; pointsVec[4][1]= xc1[1] + std::sin(M_PI_4)*rCore;
    pointsVec[5][0]= xc1[0] - std::cos(M_PI_4)*rCore; pointsVec[5][1]= xc1[1] + std::sin(M_PI_4)*rCore;
    pointsVec[6][0]= xc1[0] - std::cos(M_PI_4)*rCore; pointsVec[6][1]= xc1[1] - std::sin(M_PI_4)*rCore;
    pointsVec[7][0]= xc1[0] + std::cos(M_PI_4)*rCore; pointsVec[7][1]= xc1[1] - std::sin(M_PI_4)*rCore;

    //hole 2
    pointsVec[8][0]= xc2[0] + 0.5 * std::cos(M_PI_4)*rCore; pointsVec[8][1]= xc2[1] + 0.5 * std::sin(M_PI_4)*rCore;
    pointsVec[9][0]= xc2[0] - 0.5 * std::cos(M_PI_4)*rCore; pointsVec[9][1]= xc2[1] + 0.5 * std::sin(M_PI_4)*rCore;
    pointsVec[10][0]= xc2[0] - 0.5 * std::cos(M_PI_4)*rCore; pointsVec[10][1]= xc2[1];
    pointsVec[11][0]= xc2[0] + 0.5 * std::cos(M_PI_4)*rCore; pointsVec[11][1]= xc2[1];

    pointsVec[12][0]= xc2[0] + std::cos(M_PI_4)*rCore; pointsVec[12][1]= xc2[1] + std::sin(M_PI_4)*rCore;
    pointsVec[13][0]= xc2[0] - std::cos(M_PI_4)*rCore; pointsVec[13][1]= xc2[1] + std::sin(M_PI_4)*rCore;
    pointsVec[14][0]= xc2[0] - rCore; pointsVec[14][1]= xc2[1];
    pointsVec[15][0]= xc2[0] + rCore; pointsVec[15][1]= xc2[1];

    //cladding
    pointsVec[16][0]= 0.; pointsVec[16][1]= 0.;
    pointsVec[17][0]= pointsVec[6][0]; pointsVec[17][1]= 0.;
    pointsVec[18][0]= bound; pointsVec[18][1]= 0.;

    pointsVec[19][0]= pointsVec[16][0]; pointsVec[19][1]= pointsVec[12][1];
    pointsVec[20][0]= pointsVec[17][0]; pointsVec[20][1]= pointsVec[12][1];
    pointsVec[21][0]= pointsVec[18][0]; pointsVec[21][1]= pointsVec[12][1];

    pointsVec[22][0]= pointsVec[16][0]; pointsVec[22][1]= pointsVec[6][1];
    pointsVec[23][0]= pointsVec[12][0]; pointsVec[23][1]= pointsVec[6][1];
    pointsVec[24][0]= pointsVec[18][0]; pointsVec[24][1]= pointsVec[6][1];

    pointsVec[25][0]= pointsVec[16][0]; pointsVec[25][1]= pointsVec[4][1];
    pointsVec[26][0]= pointsVec[12][0]; pointsVec[26][1]= pointsVec[4][1];
    pointsVec[27][0]= pointsVec[18][0]; pointsVec[27][1]= pointsVec[4][1];

    pointsVec[28][0]= pointsVec[16][0]; pointsVec[28][1]= bound;
    pointsVec[29][0]= pointsVec[5][0]; pointsVec[29][1]= bound;
    pointsVec[30][0]= pointsVec[4][0]; pointsVec[30][1]= bound;
    pointsVec[31][0]= pointsVec[12][0]; pointsVec[31][1]= bound;
    pointsVec[32][0]= pointsVec[18][0]; pointsVec[32][1]= bound;
    //pml yp
    pointsVec[33][0]= pointsVec[16][0]; pointsVec[33][1]= bound + lengthPML;
    pointsVec[34][0]= pointsVec[5][0]; pointsVec[34][1]= bound + lengthPML;
    pointsVec[35][0]= pointsVec[4][0]; pointsVec[35][1]= bound + lengthPML;
    pointsVec[36][0]= pointsVec[12][0]; pointsVec[36][1]= bound + lengthPML;
    pointsVec[37][0]= pointsVec[18][0]; pointsVec[37][1]= bound + lengthPML;
    //pml xpyp
    pointsVec[38][0]= bound + lengthPML; pointsVec[38][1]= bound + lengthPML;
    //pml xp
    pointsVec[39][0]= bound + lengthPML; pointsVec[39][1]= pointsVec[28][1];
    pointsVec[40][0]= bound + lengthPML; pointsVec[40][1]= pointsVec[4][1];
    pointsVec[41][0]= bound + lengthPML; pointsVec[41][1]= pointsVec[6][1];
    pointsVec[42][0]= bound + lengthPML; pointsVec[42][1]= pointsVec[12][1];
    pointsVec[43][0]= bound + lengthPML; pointsVec[43][1]= pointsVec[14][1];

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
    quadPointsVec(18,0) = 27; quadPointsVec(18,1) = 26; quadPointsVec(18,2) = 23; quadPointsVec(18,3) = 24;
    quadPointsVec(19,0) = 29; quadPointsVec(19,1) = 28; quadPointsVec(19,2) = 25; quadPointsVec(19,3) = 5;
    quadPointsVec(20,0) = 30; quadPointsVec(20,1) = 29; quadPointsVec(20,2) = 5; quadPointsVec(20,3) = 4;
    quadPointsVec(21,0) = 31; quadPointsVec(21,1) = 30; quadPointsVec(21,2) = 4; quadPointsVec(21,3) = 26;
    quadPointsVec(22,0) = 32; quadPointsVec(22,1) = 31; quadPointsVec(22,2) = 26; quadPointsVec(22,3) = 27;
    ///////////////////PML XP
    quadPointsVec(23,0) = 39; quadPointsVec(23,1) = 32; quadPointsVec(23,2) = 27; quadPointsVec(23,3) = 40;
    quadPointsVec(24,0) = 40; quadPointsVec(24,1) = 27; quadPointsVec(24,2) = 24; quadPointsVec(24,3) = 41;
    quadPointsVec(25,0) = 41; quadPointsVec(25,1) = 24; quadPointsVec(25,2) = 21; quadPointsVec(25,3) = 42;
    quadPointsVec(26,0) = 42; quadPointsVec(26,1) = 21; quadPointsVec(26,2) = 18; quadPointsVec(26,3) = 43;
    ///////////////////PML YP
    quadPointsVec(27,0) = 34; quadPointsVec(27,1) = 33; quadPointsVec(27,2) = 28; quadPointsVec(27,3) = 29;
    quadPointsVec(28,0) = 35; quadPointsVec(28,1) = 34; quadPointsVec(28,2) = 29; quadPointsVec(28,3) = 30;
    quadPointsVec(29,0) = 36; quadPointsVec(29,1) = 35; quadPointsVec(29,2) = 30; quadPointsVec(29,3) = 31;
    quadPointsVec(30,0) = 37; quadPointsVec(30,1) = 36; quadPointsVec(30,2) = 31; quadPointsVec(30,3) = 32;
    ///////////////////PML XPYP
    quadPointsVec(31,0) = 38; quadPointsVec(31,1) = 37; quadPointsVec(31,2) = 32; quadPointsVec(31,3) = 39;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    // quadrilateral zero to four - first hole                                                              //
    // quadrilateral five to eight - second hole                                                            //
    // quadrilateral nine to twenty-eight - fiber's cladding                                                //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////CREATE QUADS///////////////////////////////////////
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
    matIdsQuads[23] = matIdPMLxp;    matIdsQuads[24] = matIdPMLxp;
    matIdsQuads[25] = matIdPMLxp;    matIdsQuads[26] = matIdPMLxp;
    matIdsQuads[27] = matIdPMLyp;    matIdsQuads[28] = matIdPMLyp;
    matIdsQuads[29] = matIdPMLyp;    matIdsQuads[30] = matIdPMLyp;
    matIdsQuads[31] = matIdPMLxpyp;

    TPZManVector<int,nQuads> nDivQsi(nQuads,-1);
    TPZManVector<int,nQuads> nDivEta(nQuads,-1);
    //first hole
    nDivQsi[0]=nDivRHole;nDivEta[0]=nDivRHole;
    nDivQsi[1]=nDivRHole; nDivEta[1]=nDivTHole;
    nDivQsi[2]=nDivTHole; nDivEta[2]=nDivRHole;
    nDivQsi[3]=nDivTHole; nDivEta[3]=nDivRHole;
    nDivQsi[4]=nDivRHole; nDivEta[4]=nDivTHole;
    //second hole
    nDivQsi[5]=nDivRHole; nDivEta[5]=(int)std::ceil(nDivRHole/2.);
    nDivQsi[6]=nDivRHole; nDivEta[6]=nDivTHole;
    nDivQsi[7]=nDivTHole; nDivEta[7]=(int)std::ceil(nDivRHole/2.);
    nDivQsi[8]=nDivTHole; nDivEta[8]=(int)std::ceil(nDivRHole/2.);
    //cladding
    nDivQsi[9]= nDivCladding3; nDivEta[9]=(int)std::ceil(nDivRHole/2.);
    nDivQsi[10]=nDivRHole; nDivEta[10]=(int)std::ceil(nDivRHole/2.);
    nDivQsi[11]=nDivCladding1; nDivEta[11]=(int)std::ceil(nDivRHole/2.);

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

    nDivQsi[23]=nDivPml; nDivEta[23]=nDivCladding1;
    nDivQsi[24]=nDivPml; nDivEta[24]=nDivRHole;
    nDivQsi[25]=nDivPml; nDivEta[25]=nDivCladding2;
    nDivQsi[26]=nDivPml; nDivEta[26]=(int)std::ceil(nDivRHole/2.);

    nDivQsi[27]=nDivCladding3; nDivEta[27]=nDivPml;
    nDivQsi[28]=nDivRHole; nDivEta[28]=nDivPml;
    nDivQsi[29]=nDivRHole; nDivEta[29]=nDivPml;
    nDivQsi[30]=nDivCladding1; nDivEta[30]=nDivPml;

    nDivQsi[31]=nDivPml; nDivEta[31]=nDivPml;

    TPZManVector<bool,nQuads> side1NonLinearVec(nQuads,false);
    TPZManVector<bool,nQuads> side2NonLinearVec(nQuads,false);
    TPZManVector<bool,nQuads> side3NonLinearVec(nQuads,false);
    TPZManVector<bool,nQuads> side4NonLinearVec(nQuads,false);

    side1NonLinearVec[1] = true;
    side2NonLinearVec[2] = true;
    side4NonLinearVec[3] = true;
    side3NonLinearVec[4] = true;

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


    TPZManVector<TPZVec<REAL>,14> thetaVec(14,TPZVec<REAL>(2,0.));

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

    TPZVec<int> matIdBoundVec(4,-1);
    matIdBoundVec[0] = matIdSymmetryX;//low
    matIdBoundVec[1] = matIdBC;//down
    matIdBoundVec[2] = matIdBC;//up
    matIdBoundVec[3] = matIdSymmetryY;//left

    TPZVec<REAL> boundDistVec(4,-1);
    boundDistVec[0] = 0.;//low
    boundDistVec[1] = bound+lengthPML;//down
    boundDistVec[2] = bound+lengthPML;//up
    boundDistVec[3] = 0.;//left

    TPZVec<REAL> rVec(14,rCore);

    gmesh = CreateStructuredMesh(pointsVec, edgesVec, quadPointsVec, matIdsQuads, nDivQsi,
                         nDivEta, side1NonLinearVec, side2NonLinearVec, side3NonLinearVec, side4NonLinearVec,
                         thetaVec, xcRef, matIdBoundVec, boundDistVec, rVec);

    matIdVec.Resize(8);
    matIdVec[0]=matIdHole;
    matIdVec[1]=matIdCladding;
    pmlTypeVec.Resize(3);
    pmlTypeVec[0]=SPZModalAnalysisData::xp;
    matIdVec[2]=matIdPMLxp;
    pmlTypeVec[1]=SPZModalAnalysisData::yp;
    matIdVec[3]=matIdPMLyp;
    pmlTypeVec[2]=SPZModalAnalysisData::xpyp;
    matIdVec[4]=matIdPMLxpyp;

    boundTypeVec.Resize(3);
    boundTypeVec[0] = SPZModalAnalysisData::PEC;
    matIdVec[matIdVec.size()-3]=matIdBC;
    boundTypeVec[1] = symmetryX;
    matIdVec[matIdVec.size()-2]=matIdSymmetryX;
    boundTypeVec[2] = symmetryY;
    matIdVec[matIdVec.size()-1]=matIdSymmetryY;

    gmesh->BuildConnectivity();

    TPZVec<TPZGeoEl *> sons;


    if(refine){
        gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
        auto uniformref = gRefDBase.GetUniformRefPattern(EQuadrilateral);
        int nel = gmesh->NElements();
        for(int i = 0; i< nel; i++){
            TPZGeoEl * gel = gmesh->Element(i);
            TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad>> *geoBlend =
                    dynamic_cast<TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad>> *> (gmesh->Element(i));
            bool condMat = (gel->MaterialId() == matIdCladding);
            if(gel->Type() == EOned) continue;
            TPZVec<REAL> qsiPos(2,0.);
            TPZVec<REAL> xPos;
            gel->X(qsiPos,xPos);
            int nRefines = 0;
            for(int iRule = 0; iRule < refineRules.size(); iRule++){
                nRefines = refineRules[iRule](xPos) ? nRefines+1 : nRefines;
            }
            auto oldref = gel->GetRefPattern();
            TPZVec<TPZGeoEl *> thisTimeRefined(1,gel);
            for(int refLevel = 0; refLevel < nRefines; refLevel ++){
                TPZVec<TPZGeoEl *> nextTimeRefined(0,nullptr);
                for(int iEl = 0; iEl < thisTimeRefined.size(); iEl++){
                    TPZGeoEl *geo = thisTimeRefined[iEl];
                    geo->X(qsiPos,xPos);
                    int shouldRefine = 0;
                    for(int iref = 0; iref < refineRules.size(); iref++){
                        shouldRefine = refineRules[iref](xPos) ? shouldRefine+1 : shouldRefine;
                    }
                    if(shouldRefine < refLevel) continue;
                    auto oldref = geo->GetRefPattern();
                    geo->SetRefPattern(uniformref);
                    geo->Divide(sons);
                    const int oldSize = nextTimeRefined.size();
                    nextTimeRefined.Resize(oldSize+sons.size());
                    for(int i = 0; i < sons.size(); i++){
                        sons[i]->SetRefPattern(oldref);
                        nextTimeRefined[oldSize+i] = sons[i];
                    }
                }
                thisTimeRefined = nextTimeRefined;
            }
        }
        gmesh->ResetConnectivities();
        gmesh->BuildConnectivity();
    }

    long nel = gmesh->NElements();
    for (long iel = 0; iel < nel; iel++) {
        TPZGeoEl *gelQuad = gmesh->ElementVec()[iel];
        if(gelQuad->Type() == EQuadrilateral && gelQuad->HasSubElement() == 0){
            gelQuad->Divide(sons); nel = gmesh->NElements(); continue;
        }
    }
    gmesh->ResetConnectivities();
    gmesh->BuildConnectivity();


    if(print){
        std::string meshFileName = prefix + "gmesh";
        const size_t strlen = meshFileName.length();
        meshFileName.append(".vtk");
        std::ofstream outVTK(meshFileName.c_str());
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outVTK, true);
        outVTK.close();
        meshFileName.replace(strlen, 4, ".txt");
        std::ofstream outTXT(meshFileName.c_str());
        gmesh->Print(outTXT);
        outTXT.close();
    }

//    #ifdef PZDEBUG
//    TPZCheckGeom * Geometrytest = new TPZCheckGeom(gmesh);
//    int isBadMeshQ = Geometrytest->PerformCheck();
//
////	if (isBadMeshQ) {
////		DebugStop();
////	}
//    #endif
//    return;
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

void CreateCMesh(TPZVec<TPZCompMesh *> &meshVecOut, TPZGeoMesh *gmesh, int pOrder, const TPZVec<int> &matIdVec,
                 TPZVec<STATE> &urVec, TPZVec<STATE> &erVec, REAL f0, bool isCutOff, const std::string &prefix,
                 const bool &print, const REAL &scale, const REAL &alphaMax,
                 TPZVec<SPZModalAnalysisData::boundtype> &boundTypeVec,
                 TPZVec<SPZModalAnalysisData::pmltype> &pmlTypeVec, bool refineP,
                 TPZVec<std::function<bool (const TPZVec<REAL> &)>> refineRules, SPZModalAnalysisData::NedEl elType) {
    bool usingNedelecTypeTwo = elType == SPZModalAnalysisData::NedEl::TypeTwo ? true : false;
    const int dim = 2;

    TPZManVector<int, 8> volMatIdVec(matIdVec.size() - boundTypeVec.size() - pmlTypeVec.size());
    for (int i = 0; i < volMatIdVec.size(); i++) {
        volMatIdVec[i] = matIdVec[i];
    }
    TPZManVector<int, 8> pmlMatIdVec(pmlTypeVec.size());
    for (int i = 0; i < pmlMatIdVec.size(); i++) {
        pmlMatIdVec[i] = matIdVec[volMatIdVec.size() + i];
    }
    TPZManVector<int, 8> boundMatIdVec(boundTypeVec.size());
    for (int i = 0; i < boundMatIdVec.size(); i++) {
        boundMatIdVec[i] = matIdVec[volMatIdVec.size() + pmlMatIdVec.size() + i];
    }

    if (volMatIdVec.size() != urVec.size()) {
        std::cout << "The number of materials is not consistent with the mesh materials." << std::endl;
        std::cout << "The number of materials is " << urVec.size() << " and there are " << volMatIdVec.size();
        std::cout << " mesh materials." << std::endl;
        DebugStop();
    }
    if (pmlMatIdVec.size() != pmlTypeVec.size()) {
        std::cout << "The number of PML types is not consistent with the mesh PMLs." << std::endl;
        std::cout << "The number of PML types is " << pmlTypeVec.size() << " and there are " << pmlMatIdVec.size();
        std::cout << " PMLs in the mesh." << std::endl;
        DebugStop();
    }
    if (boundMatIdVec.size() != boundTypeVec.size()) {
        std::cout << "The number of Boundary Conditions types is not consistent with the mesh." << std::endl;
        std::cout << "The number of BC types is " << boundTypeVec.size() << " and there are " << boundMatIdVec.size();
        std::cout << " boundaries in the mesh." << std::endl;
        DebugStop();
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

    TPZCompMesh *cmeshH1 = new TPZCompMesh(gmesh);
    if(usingNedelecTypeTwo) cmeshH1->SetDefaultOrder(pOrder + 1); // seta ordem polimonial de aproximacao
    else cmeshH1->SetDefaultOrder(pOrder); // seta ordem polimonial de aproximacao
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
    for (int i = 0; i < boundMatIdVec.size(); i++) {
        bcond = matH1->CreateBC(matH1, boundMatIdVec[i], boundTypeVec[i], val1, val2);
        cmeshH1->InsertMaterialObject(bcond);
    }

    cmeshH1->SetAllCreateFunctionsContinuous();
    cmeshH1->AutoBuild(set);
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

    for (int i = 0; i < boundMatIdVec.size(); i++) {
        bcond = matHCurl->CreateBC(matHCurl, boundMatIdVec[i], boundTypeVec[i], val1, val2);
        cmeshHCurl->InsertMaterialObject(bcond);
    }

    if(usingNedelecTypeTwo) {
        cmeshHCurl->SetAllCreateFunctionsHDiv();
        cmeshHCurl->AutoBuild(set);
        cmeshHCurl->CleanUpUnconnectedNodes();
        int64_t nel = cmeshHCurl->NElements();
        for (int64_t el = 0; el<nel; el++) {
            TPZCompEl *cel = cmeshHCurl->Element(el);
            if(!cel) continue;
            TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
            int nc = intel->NConnects();
            if (nc <=1) {
                continue;
            }
            TPZGeoEl *gel = intel->Reference();
            int ns = gel->NSides();
            intel->ForceSideOrder(ns-1, pOrder -1);
        }
        cmeshHCurl->ExpandSolution();
    }
    else {
        cmeshHCurl->SetAllCreateFunctionsHCurl(); // define espaco de aproximacao
        cmeshHCurl->AutoBuild(set);
        cmeshHCurl->CleanUpUnconnectedNodes();
    }

    TPZCompMesh *cmeshMF = new TPZCompMesh(gmesh);
    TPZManVector<TPZMatModalAnalysis *, 8> matMultiPhysics(volMatIdVec.size() + pmlMatIdVec.size());

    if (isCutOff) {
        for (int i = 0; i < volMatIdVec.size(); ++i) {
            matMultiPhysics[i] = new TPZMatWaveguideCutOffAnalysis(volMatIdVec[i], f0, urVec[i], erVec[i], 1. / scale);
            cmeshMF->InsertMaterialObject(matMultiPhysics[i]);
        }
    } else {
        for (int i = 0; i < volMatIdVec.size(); ++i) {
            if(usingNedelecTypeTwo) matMultiPhysics[i] =
                    new TPZMatModalAnalysisHDiv(volMatIdVec[i], f0, urVec[i], erVec[i], 1. / scale);
            else matMultiPhysics[i] = new TPZMatModalAnalysis(volMatIdVec[i], f0, urVec[i], erVec[i], 1. / scale);
            cmeshMF->InsertMaterialObject(matMultiPhysics[i]);
        }
    }
    for (int i = 0; i < pmlMatIdVec.size(); i++) {
        REAL xMax = -1e20, xMin = 1e20, yMax = -1e20, yMin = 1e20;
        TPZGeoMesh *gmesh = cmeshMF->Reference();
        for (int iel = 0; iel < gmesh->NElements(); ++iel) {
            TPZGeoEl *geo = gmesh->Element(iel);
            if (geo->MaterialId() == pmlMatIdVec[i]) {
                for (int iNode = 0; iNode < geo->NCornerNodes(); ++iNode) {
                    TPZManVector<REAL, 3> co(3);
                    geo->Node(iNode).GetCoordinates(co);
                    const REAL &xP = co[0];
                    const REAL &yP = co[1];
                    if (xP > xMax) {
                        xMax = xP;
                    }
                    if (xP < xMin) {
                        xMin = xP;
                    }
                    if (yP > yMax) {
                        yMax = yP;
                    }
                    if (yP < yMin) {
                        yMin = yP;
                    }
                }
            }
        }
        bool attx, atty;
        REAL xBegin, yBegin, d;
        REAL boundPosX = -666, boundPosY = -666;
        switch (pmlTypeVec[i]) {
            case SPZModalAnalysisData::xp:
                attx = true;
                atty = false;
                xBegin = xMin;
                yBegin = -1;
                d = xMax - xMin;
                boundPosX = xBegin;
                boundPosY = (yMax + yMin)/2;
                break;
            case SPZModalAnalysisData::yp:
                attx = false;
                atty = true;
                xBegin = -1;
                yBegin = yMin;
                d = yMax - yMin;
                boundPosX = (xMax + xMin)/2;
                boundPosY = yBegin;
                break;
            case SPZModalAnalysisData::xm:
                attx = true;
                atty = false;
                xBegin = xMax;
                yBegin = -1;
                d = xMax - xMin;
                boundPosX = xBegin;
                boundPosY = (yMax + yMin)/2;
                break;
            case SPZModalAnalysisData::ym:
                attx = false;
                atty = true;
                xBegin = -1;
                yBegin = yMax;
                d = yMax - yMin;
                boundPosX = (xMax + xMin)/2;
                boundPosY = yBegin;
                break;
            case SPZModalAnalysisData::xpyp:
                attx = true;
                atty = true;
                xBegin = xMin;
                yBegin = yMin;
                d = xMax - xMin;
                boundPosX = xBegin;
                boundPosY = yBegin;
                break;
            case SPZModalAnalysisData::xmyp:
                attx = true;
                atty = true;
                xBegin = xMax;
                yBegin = yMin;
                d = xMax - xMin;
                boundPosX = xBegin;
                boundPosY = yBegin;
                break;
            case SPZModalAnalysisData::xmym:
                attx = true;
                atty = true;
                xBegin = xMax;
                yBegin = yMax;
                d = xMax - xMin;
                boundPosX = xBegin;
                boundPosY = yBegin;
                break;
            case SPZModalAnalysisData::xpym:
                attx = true;
                atty = true;
                xBegin = xMin;
                yBegin = yMax;
                d = xMax - xMin;
                boundPosX = xBegin;
                boundPosY = yBegin;
                break;
        }
        TPZGeoEl * closestEl = nullptr;
        REAL dist = 1e16;
        for(int iEl = 0; iEl < gmesh->NElements(); iEl++){
            TPZGeoEl * currentEl = gmesh->ElementVec()[iEl];
            if ( currentEl->NSubElements() > 0  || ( currentEl->Dimension() != 2 ) ) continue;
            bool isPmlMat = false;
            for(int iPml = 0; iPml < pmlMatIdVec.size(); iPml++){
                if ( currentEl->MaterialId() == pmlMatIdVec[iPml] ){
                    isPmlMat = true;
                }
            }
            if(isPmlMat) continue;
            TPZVec<REAL> qsi(2,-1);
            const int largerSize = currentEl->NSides() - 1;
            currentEl->CenterPoint(largerSize, qsi);
            TPZVec<REAL> xCenter(3,-1);
            currentEl->X(qsi, xCenter);
            const REAL currentDist = (xCenter[0]-boundPosX)*(xCenter[0]-boundPosX) +
                    (xCenter[1]-boundPosY)*(xCenter[1]-boundPosY);
            if(currentDist < dist){
                dist = currentDist;
                closestEl = currentEl;
            }
        }
        const int pmlNeighbourMatId = closestEl->MaterialId();
        int outerMaterialPos = -1;
        bool matFound = false;
        for(int imat = 0; imat < volMatIdVec.size(); imat++){
            if(matMultiPhysics[imat] && matMultiPhysics[imat]->Id() == pmlNeighbourMatId){
                outerMaterialPos = imat;
                matFound = true;
                break;
            }
        }
        if(matFound == false){
            PZError<<"could not find pml neighbour, aborting...."<<std::endl;
            DebugStop();
        }
        if(usingNedelecTypeTwo){
            matMultiPhysics[volMatIdVec.size() + i] =
                    new TPZMatWaveguidePmlHDiv(pmlMatIdVec[i],
                                               *matMultiPhysics[outerMaterialPos],
                                               attx, xBegin,
                                               atty, yBegin,
                                               alphaMax, d);
        }
        else{
            matMultiPhysics[volMatIdVec.size() + i] =
                    new TPZMatWaveguidePml(pmlMatIdVec[i],
                                           *matMultiPhysics[outerMaterialPos],
                                           attx, xBegin,
                                           atty, yBegin,
                                           alphaMax, d);
        }
        cmeshMF->InsertMaterialObject(matMultiPhysics[volMatIdVec.size() + i]);
    }


    for (int i = 0; i < boundMatIdVec.size(); i++) {
        bcond = matMultiPhysics[0]->CreateBC(matMultiPhysics[0], boundMatIdVec[i], boundTypeVec[i], val1, val2);
        cmeshMF->InsertMaterialObject(bcond);
    }

    cmeshMF->SetDimModel(dim);
    cmeshMF->SetAllCreateFunctionsMultiphysicElem();

    cmeshMF->AutoBuild(set);
    cmeshMF->CleanUpUnconnectedNodes();

    TPZVec<TPZCompMesh *> meshVec(2);
    meshVec[matMultiPhysics[0]->H1Index()] = cmeshH1;
    meshVec[matMultiPhysics[0]->HCurlIndex()] = cmeshHCurl;

    if (refineP) {
        for(int iMesh = 0; iMesh < meshVec.size(); iMesh++){
            gmesh->ResetReference();
            meshVec[iMesh]->LoadReferences();
            for (int i = 0; i < meshVec[iMesh]->NElements(); i++) {
                TPZCompEl *cel = meshVec[iMesh]->Element(i);
                TPZGeoEl *gel = cel->Reference();
                if (gel->Father()) {
                    do{
                        gel = gel->Father();
                    }
                    while (gel->Father() && gel->Father()->Type() != EQuadrilateral);
                }
                TPZVec<REAL> qsiPos(2, 0.);
                TPZVec<REAL> xPos;
                gel->X(qsiPos, xPos);
                int actualPOrder = pOrder;
                for (int iRule = 0; iRule < refineRules.size(); iRule++) {
                    bool res = refineRules[iRule](xPos);
                    if (res) actualPOrder++;
                }
                TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
                actualPOrder = actualPOrder > 5 ? 5 : actualPOrder;
                if(pOrder != actualPOrder) intel->PRefine(actualPOrder);
            }
            meshVec[iMesh]->ExpandSolution();
        }
    }

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
