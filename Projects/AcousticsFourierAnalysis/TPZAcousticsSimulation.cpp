#include "TPZAcousticsSimulation.h"
#include "TPZAcousticGeoMesher.h"
#include "TPZAcousticCompMesher.h"
#include "TPZAcousticAnalysis.h"
#include "TPZAcousticFreqDomainAnalysis.h"
#include "TPZAcousticTimeDomainAnalysis.h"
#include "pzcheckgeom.h"
#include "pzcmesh.h"
#include "TPZGmshReader.h"
#include "TPZVTKGeoMesh.h"
#include "pzbndcond.h"
#include "pzgeoelrefless.h"
#include <pzgeopoint.h>
#include <boost/date_time/posix_time/posix_time.hpp>

void TPZAcousticsSimulation::RunSimulation() {
    //////////////////////////////////////////////////////
    //////////////////// CREATE GMESH ////////////////////
    //////////////////////////////////////////////////////
    boost::posix_time::ptime t1_total =
            boost::posix_time::microsec_clock::local_time();
    std::cout << "Creating gmesh... ";
    boost::posix_time::ptime t1_g =
            boost::posix_time::microsec_clock::local_time();

    std::string meshFileName = this->fSimData.fSimulationSettings.meshName;
    std::string prefix = this->fSimData.fOutputSettings.resultsDir;
    TPZAcousticGeoMesher geoMesh(meshFileName, prefix);
    {

        REAL nElemPerLambdaTimesOmega = this->fSimData.fSimulationSettings.nElemPerLambda *
                                        this->fSimData.fSourceSettings.centralFrequency;
        geoMesh.CreateGMesh(nElemPerLambdaTimesOmega);
        bool &printGtxt = this->fSimData.fOutputSettings.printGmeshTxt;
        bool &printGvtk = this->fSimData.fOutputSettings.printGmeshVtk;
        std::string outFilesName = "geomesh";
        geoMesh.PrintMesh(outFilesName, prefix, printGvtk, printGtxt);
    }
    ////////////////////////////////////////////////////////////////////////
    //////////////////////////CREATE SOURCE GEO EL//////////////////////////
    ////////////////////////////////////////////////////////////////////////


    geoMesh.CreateSourceNode(this->fSimData.fSourceSettings.posX, this->fSimData.fSourceSettings.posY);


    boost::posix_time::ptime t2_g =
            boost::posix_time::microsec_clock::local_time();
    std::cout << "Created! " << t2_g - t1_g << std::endl;

    //////////////////////////////////////////////////////
    //////////////////// CREATE CMESH ////////////////////
    //////////////////////////////////////////////////////
    std::cout << "Creating cmesh... ";
    boost::posix_time::ptime t1_c =
            boost::posix_time::microsec_clock::local_time();
    bool isAxisymmetric = this->fSimData.fSimulationSettings.axiSymmetricSimulation;
    auto boundTypeVec = this->fSimData.fSimulationSettings.boundType;
    TPZAcousticCompMesher compMesh(&geoMesh, boundTypeVec, isAxisymmetric);
    {
        const int &pOrder = this->fSimData.fSimulationSettings.pOrder;
        if(this->fSimData.fSimulationSettings.simType == SPZAcousticData::ESimulationType::frequencyDomain){
            compMesh.CreateFourierMesh(pOrder);
        }else{
            compMesh.CreateTransientMesh(pOrder);
        }
        std::string outFilesName = "compmesh";
        const bool &printCtxt = this->fSimData.fOutputSettings.printCmeshTxt;
        const bool &printCvtk = this->fSimData.fOutputSettings.printCmeshVtk;
        compMesh.PrintMesh(outFilesName,prefix,printCvtk,printCtxt);
    }
    boost::posix_time::ptime t2_c =
            boost::posix_time::microsec_clock::local_time();
    std::cout << "Created! " << t2_c - t1_c << std::endl;


    //////////////////////////////////////////////////////
    ////////////////// CALCULATE DELTA T /////////////////
    //////////////////////////////////////////////////////
    auto elSizeMap = geoMesh.GetElSizes();
    auto velocityMap = geoMesh.GetVelocityMap();
    const REAL totalTime = this->fSimData.fSimulationSettings.totalTime;
    REAL deltaT = -1;
    int nTimeSteps = -1;

    if(this->fSimData.fSimulationSettings.isCflBound){
        deltaT = CalculateDeltaT(elSizeMap, velocityMap);
        nTimeSteps = (int) std::ceil(totalTime/deltaT) + 1;
        deltaT = totalTime / (nTimeSteps - 1);
    }else{
        nTimeSteps = this->fSimData.fSimulationSettings.nTimeSteps;
        deltaT = totalTime / (nTimeSteps - 1);
    }

    const int &nThreads = this->fSimData.fSimulationSettings.nThreads;
    TPZAcousticAnalysis *analysis = nullptr;

    bool &filter = this->fSimData.fSimulationSettings.filterBoundaryEqs;
    if(this->fSimData.fSimulationSettings.simType == SPZAcousticData::ESimulationType::frequencyDomain){
        analysis = new TPZAcousticFreqDomainAnalysis(&compMesh, nThreads, deltaT, nTimeSteps, filter);
    }else{
        analysis = new TPZAcousticTimeDomainAnalysis(&compMesh, nThreads,deltaT, nTimeSteps, filter);
    }

    analysis->InitializeComputations();
    {
        const REAL &wZero = this->fSimData.fSourceSettings.centralFrequency;
        const REAL &peakTime = this->fSimData.fSourceSettings.peakTime;
        const REAL &amplitude = this->fSimData.fSourceSettings.amplitude;
        analysis->SetUpGaussianSource(wZero, peakTime, amplitude);
    }

    if(this->fSimData.fSimulationSettings.simType == SPZAcousticData::ESimulationType::frequencyDomain){
        const REAL & wMax = this->fSimData.fFourierSettings.wMax;
        const int & nSamples = this->fSimData.fFourierSettings.nSamples;
        const REAL & alphaFreqShift = this->fSimData.fFourierSettings.alphaFreqShift;
        TPZAcousticFreqDomainAnalysis *dummyAnalysis = dynamic_cast<TPZAcousticFreqDomainAnalysis *>(analysis);
        dummyAnalysis->SetUpFourierSettings(wMax,nSamples,alphaFreqShift);
    }

    analysis->RunSimulationSteps();
    boost::posix_time::ptime t2_total =
            boost::posix_time::microsec_clock::local_time();
    std::cout << "Finished computations! " << t2_total - t1_total << std::endl;
    if(this->fSimData.fOutputSettings.vtkSol)
    {
        const int vtkRes = this->fSimData.fOutputSettings.vtkResolution;

        analysis->PostProcess(vtkRes,prefix);
    }
    delete analysis;
}

//void TPZAcousticsSimulation::FilterBoundaryEquations(TPZCompMesh *cmesh, TPZVec<int64_t> &activeEquations, int &neq, int &neqOriginal){
//    TPZManVector<int64_t, 1000> allConnects;
//    std::set<int64_t> boundConnects;
//
//    for (int iel = 0; iel < cmesh->NElements(); iel++) {
//        TPZCompEl *cel = cmesh->ElementVec()[iel];
//        if (cel == NULL) {
//            continue;
//        }
//        if (cel->Reference() == NULL) {
//            continue;
//        }
//        TPZBndCond *mat = dynamic_cast<TPZBndCond *>(cmesh->MaterialVec()[cel->Reference()->MaterialId()]);
//        if (mat && mat->Type() == 0) {
//            std::set<int64_t> boundConnectsEl;
//            cel->BuildConnectList(boundConnectsEl);
//
//            for (std::set<int64_t>::iterator iT = boundConnectsEl.begin();
//                 iT != boundConnectsEl.end(); iT++) {
//                const int64_t val = *iT;
//                if (boundConnects.find(val) == boundConnects.end()) {
//                    boundConnects.insert(val);
//                }
//            }
//        }
//    }
//    neqOriginal = cmesh->NEquations();
//    for (int iCon = 0; iCon < cmesh->NConnects(); iCon++) {
//        if (boundConnects.find(iCon) == boundConnects.end()) {
//            TPZConnect &con = cmesh->ConnectVec()[iCon];
//            int seqnum = con.SequenceNumber();
//            int pos = cmesh->Block().Position(seqnum);
//            int blocksize = cmesh->Block().Size(seqnum);
//            if (blocksize == 0)
//                continue;
//
//            int vs = activeEquations.size();
//            activeEquations.Resize(vs + blocksize);
//            for (int ieq = 0; ieq < blocksize; ieq++) {
//                activeEquations[vs + ieq] = pos + ieq;
//                neq++;
//            }
//        }
//    }
//    neqOriginal = cmesh->NEquations();
//    std::cout << "# equations(before): " << neqOriginal << std::endl;
//    std::cout << "# equations(after): " << neq << std::endl;
//    return;
//}

REAL TPZAcousticsSimulation::CalculateDeltaT(std::map<int,REAL> elSizeMap, std::map<int,REAL> velocityMap) {
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
        REAL smallerTimeStep = 1e12;
        for(auto iVelocity : velocityMap){
            auto matId = iVelocity.first;
            elSizeMap[matId] = cflVal * elSizeMap[matId] / iVelocity.second;
            if(elSizeMap[matId] < smallerTimeStep) smallerTimeStep = elSizeMap[matId];
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
        for(auto iVelocity : velocityMap){
            auto matId = iVelocity.first;
            const REAL nElemPerLambdaTimesOmega = fSimData.fSimulationSettings.nElemPerLambda *
                                                  fSimData.fSourceSettings.centralFrequency;
#ifdef PZDEBUG
            std::cout<<"material "<<matId<<"\telradius "<<elSizeMap[matId]<<"\telradius_calc";
            std::cout<<2*M_PI*iVelocity.second/nElemPerLambdaTimesOmega<<"\t velocity "<<iVelocity.second<<std::endl;
#endif
            elSizeMap[matId] = iVelocity.second * (totalTime/nTimeSteps)/ elSizeMap[matId];
#ifdef PZDEBUG
            std::cout<<"CFL for material "<<matId<<" is:"<<elSizeMap[matId]<<std::endl;
#endif
        }
    }
    const REAL deltaT = totalTime/nTimeSteps;
    std::cout<<"delta t = "<<deltaT<<std::endl;
    return deltaT;
}
