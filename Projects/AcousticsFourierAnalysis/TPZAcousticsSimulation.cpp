#include "TPZAcousticsSimulation.h"
#include "TPZAcousticGeoMesher.h"
#include "pzcheckgeom.h"
#include "pzcmesh.h"
#include "TPZGmshReader.h"
#include "TPZVTKGeoMesh.h"
#include "pzbndcond.h"
#include "pzgeoelrefless.h"
#include <pzgeopoint.h>
#include <boost/date_time/posix_time/posix_time.hpp>

void TPZAcousticsSimulation::RunSimulation(){
    //////////////////////////////////////////////////////
    //////////////////// CREATE GMESH ////////////////////
    //////////////////////////////////////////////////////
    std::cout<<"Creating gmesh... ";
    boost::posix_time::ptime t1_g =
            boost::posix_time::microsec_clock::local_time();
    TPZGeoMesh *gmesh = nullptr;
    TPZVec<int> matIdVec;
    std::map<int,REAL> rhoMap;
    std::map<int,REAL> velocityMap;
    TPZVec<REAL> elSizeVec;
    {
        std::string meshFileName = this->fSimData.fOutputSettings.resultsDir;
        std::string prefix = this->fSimData.fSimulationSettings.meshName;
        REAL nElemPerLambdaTimesOmega = this->fSimData.fSimulationSettings.nElemPerLambda *
                                        this->fSimData.fSourceSettings.centralFrequency;
        TPZAcousticGeoMesher::CreateGMesh(gmesh, meshFileName, prefix, nElemPerLambdaTimesOmega, matIdVec, elSizeVec,
                velocityMap, rhoMap);
        bool &printGtxt = this->fSimData.fOutputSettings.printGmeshTxt;
        bool &printGvtk = this->fSimData.fOutputSettings.printGmeshVtk;
        std::string outFilesName = "geomesh";
        TPZAcousticGeoMesher::PrintMesh(gmesh, outFilesName, prefix, printGvtk, printGtxt);
    }
    ////////////////////////////////////////////////////////////////////////
    //////////////////////////CREATE SOURCE GEO EL//////////////////////////
    ////////////////////////////////////////////////////////////////////////
    int matIdSource = 0;
    for (int iMat = 0; iMat< matIdVec.size(); iMat++){
        matIdSource += matIdVec[iMat];
    }
    TPZVec<int> matIdVecCp(matIdVec);
    matIdVec.resize(matIdVec.size()+1);
    matIdVec[0] = matIdSource;
    for (int iMat = 0; iMat< matIdVecCp.size(); iMat++){
        matIdVec[iMat+1] = matIdVecCp[iMat];
    }

    TPZAcousticGeoMesher::CreateSourceNode(gmesh, matIdSource, this->fSimData.fSourceSettings.posX,
            this->fSimData.fSourceSettings.posY);


    boost::posix_time::ptime t2_g =
            boost::posix_time::microsec_clock::local_time();
    std::cout<<"Created! "<<t2_g-t1_g<<std::endl;
    //////////////////////////////////////////////////////
    ////////////////// CALCULATE DELTA T /////////////////
    //////////////////////////////////////////////////////
    const REAL deltaT = CalculateDeltaT(elSizeVec, velocityMap);
    //////////////////////////////////////////////////////
    //////////////////// CREATE CMESH ////////////////////
    //////////////////////////////////////////////////////
    std::cout<<"Creating cmesh... ";
    boost::posix_time::ptime t1_c =
            boost::posix_time::microsec_clock::local_time();
    TPZCompMesh *cmesh = NULL;
    {
        const int &pOrder = this->fSimData.fSimulationSettings.pOrder;
        std::string &prefix =this->fSimData.fOutputSettings.resultsDir;
        const bool &printCtxt = this->fSimData.fOutputSettings.printCmeshTxt;
        const bool &printCvtk = this->fSimData.fOutputSettings.printCmeshVtk;

        CreateCMesh(cmesh, gmesh, pOrder, prefix, printCvtk, matIdVec, rhoMap, velocityMap,
                this->fSimData.fSimulationSettings.boundType);
    }
    boost::posix_time::ptime t2_c =
            boost::posix_time::microsec_clock::local_time();
    std::cout<<"Created! "<<t2_c-t1_c<<std::endl;


}

void TPZAcousticsSimulation::FilterBoundaryEquations(TPZCompMesh *cmesh, TPZVec<int64_t> &activeEquations, int &neq, int &neqOriginal){
    TPZManVector<int64_t, 1000> allConnects;
    std::set<int64_t> boundConnects;

    for (int iel = 0; iel < cmesh->NElements(); iel++) {
        TPZCompEl *cel = cmesh->ElementVec()[iel];
        if (cel == NULL) {
            continue;
        }
        if (cel->Reference() == NULL) {
            continue;
        }
        TPZBndCond *mat = dynamic_cast<TPZBndCond *>(cmesh->MaterialVec()[cel->Reference()->MaterialId()]);
        if (mat && mat->Type() == 0) {
            std::set<int64_t> boundConnectsEl;
            cel->BuildConnectList(boundConnectsEl);

            for (std::set<int64_t>::iterator iT = boundConnectsEl.begin();
                 iT != boundConnectsEl.end(); iT++) {
                const int64_t val = *iT;
                if (boundConnects.find(val) == boundConnects.end()) {
                    boundConnects.insert(val);
                }
            }
        }
    }
    neqOriginal = cmesh->NEquations();
    for (int iCon = 0; iCon < cmesh->NConnects(); iCon++) {
        if (boundConnects.find(iCon) == boundConnects.end()) {
            TPZConnect &con = cmesh->ConnectVec()[iCon];
            int seqnum = con.SequenceNumber();
            int pos = cmesh->Block().Position(seqnum);
            int blocksize = cmesh->Block().Size(seqnum);
            if (blocksize == 0)
                continue;

            int vs = activeEquations.size();
            activeEquations.Resize(vs + blocksize);
            for (int ieq = 0; ieq < blocksize; ieq++) {
                activeEquations[vs + ieq] = pos + ieq;
                neq++;
            }
        }
    }
    neqOriginal = cmesh->NEquations();
    std::cout << "# equations(before): " << neqOriginal << std::endl;
    std::cout << "# equations(after): " << neq << std::endl;
    return;
}

REAL TPZAcousticsSimulation::CalculateDeltaT(const TPZVec<REAL> &elSizeVec, std::map<int,REAL> velocityMap) {
    REAL &totalTime = this->fSimData.fSimulationSettings.totalTime;
    int &nTimeSteps = this->fSimData.fSimulationSettings.nTimeSteps;
    const REAL &cflVal = this->fSimData.fSimulationSettings.cfl;

    if(this->fSimData.fSimulationSettings.isCflBound)
    {
        if(cflVal< 0){
            std::cout<<"You have not set the CFL"<<std::endl;
            DebugStop();
            exit(1);
        }
        std::cout<<"Calculating time step for CFL = "<<std::setprecision(4)<<cflVal<<std::endl;
        boost::posix_time::ptime t1_c =
                boost::posix_time::microsec_clock::local_time();
        std::map<int,REAL> sizeMap;//matId, minElSize
        int i = 0;
        for(auto iVelocity : velocityMap){
            auto matId = iVelocity.first;
            sizeMap[matId] = elSizeVec[i];
            i++;
        }
        REAL smallerTimeStep = 1e12;
        for(auto iVelocity : velocityMap){
            auto matId = iVelocity.first;
            sizeMap[matId] = cflVal * sizeMap[matId] / iVelocity.second;
            if(sizeMap[matId] < smallerTimeStep) smallerTimeStep = sizeMap[matId];
        }
        nTimeSteps = std::ceil(totalTime/smallerTimeStep);
        boost::posix_time::ptime t2_c =
                boost::posix_time::microsec_clock::local_time();
        std::cout<<"took "<<t2_c-t1_c<<std::endl;
    }else{
        if(nTimeSteps < 1){
            std::cout<<"You have not set the number of time steps"<<std::endl;
            DebugStop();
            exit(1);
        }
        std::map<int,REAL> sizeMap;//matId, minElSize
        int i = 0;
        for(auto iVelocity : velocityMap){
            auto matId = iVelocity.first;
            sizeMap[matId] = elSizeVec[i];
            i++;
        }
        for(auto iVelocity : velocityMap){
            auto matId = iVelocity.first;
            const REAL nElemPerLambdaTimesOmega = fSimData.fSimulationSettings.nElemPerLambda *
                                                  fSimData.fSourceSettings.centralFrequency;
#ifdef PZDEBUG
            std::cout<<"material "<<matId<<"\telradius "<<sizeMap[matId]<<"\telradius_calc";
            std::cout<<2*M_PI*iVelocity.second/nElemPerLambdaTimesOmega<<"\t velocity "<<iVelocity.second<<std::endl;
#endif
            sizeMap[matId] = iVelocity.second * (totalTime/nTimeSteps)/ sizeMap[matId];
#ifdef PZDEBUG
            std::cout<<"CFL for material "<<matId<<" is:"<<sizeMap[matId]<<std::endl;
#endif
        }
    }
    const REAL deltaT = totalTime/nTimeSteps;
    std::cout<<"delta t = "<<deltaT<<std::endl;
    return deltaT;
}
