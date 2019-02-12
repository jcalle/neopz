#include "TPZAcousticsSimulation.h"

enum ESimulationCases{ freqAxiSymmetricHomo,freqAxiSymmetricHetero, freqHomogeneous2D, freqConcentric2D,
                      timeAxiSymmetricHomo,timeAxiSymmetricHetero, timeHomogeneous2D, timeConcentric2D};

void ConfigureConcentricCase(TPZAcousticsSimulation &sim);
void ConfigureFreqConcentricCase(TPZAcousticsSimulation &sim);

void ConfigureHomogeneousCase(TPZAcousticsSimulation &sim);
void ConfigureFreqHomogeneousCase(TPZAcousticsSimulation &sim);

void ConfigureAxiHomogeneousCase(TPZAcousticsSimulation &sim);
void ConfigureFreqAxiHomogeneousCase(TPZAcousticsSimulation &sim);

int main(int argc, char *argv[]) {
    TPZAcousticsSimulation sim;
    ESimulationCases simCase = ESimulationCases::freqAxiSymmetricHomo;
//    ESimulationCases simCase = ESimulationCases::timeAxiSymmetricHomo;
//    ESimulationCases simCase = ESimulationCases::timeHomogeneous2D;
//    ESimulationCases simCase = ESimulationCases::freqHomogeneous2D;
//    ESimulationCases simCase = ESimulationCases::freqConcentric2D;
//    ESimulationCases simCase = ESimulationCases::timeConcentric2D;
    switch(simCase){
        case ESimulationCases::freqAxiSymmetricHomo:
            ConfigureFreqAxiHomogeneousCase(sim);
            break;
        case ESimulationCases::freqHomogeneous2D:
            ConfigureFreqHomogeneousCase(sim);
            break;
        case ESimulationCases::freqConcentric2D:
            ConfigureFreqConcentricCase(sim);
            break;
        case ESimulationCases::timeAxiSymmetricHomo:
            ConfigureAxiHomogeneousCase(sim);
            break;
        case ESimulationCases::timeHomogeneous2D:
            ConfigureHomogeneousCase(sim);
            break;
        case ESimulationCases::timeConcentric2D:
            ConfigureConcentricCase(sim);
            break;
        default:
            PZError<<"What are you trying to run?"<<std::endl;
            DebugStop();
    }
    sim.RunSimulation();
}


void ConfigureConcentricCase(TPZAcousticsSimulation &sim){

    sim.fSimData.fSourceSettings.posX = 0;
    sim.fSimData.fSourceSettings.posY = 0;
    sim.fSimData.fSourceSettings.amplitude = 1e-5;
    sim.fSimData.fSourceSettings.centralFrequency = 2*M_PI*18000;
    sim.fSimData.fSourceSettings.peakTime = 1./(18000);

    sim.fSimData.fSimulationSettings.simType = SPZAcousticData::ESimulationType::timeDomain;
    sim.fSimData.fSimulationSettings.meshName = "realWEll.geo";
    sim.fSimData.fSimulationSettings.boundType.Resize(1);
    sim.fSimData.fSimulationSettings.boundType[0] = SPZAcousticData::EBoundType::softwall;
    sim.fSimData.fSimulationSettings.nElemPerLambda = 12;
    sim.fSimData.fSimulationSettings.totalTime = sim.fSimData.fSourceSettings.peakTime * 12;
    sim.fSimData.fSimulationSettings.cfl = 0.2;
    sim.fSimData.fSimulationSettings.nTimeSteps = 600;
    sim.fSimData.fSimulationSettings.isCflBound = true;
    sim.fSimData.fSimulationSettings.pOrder = 2;
    sim.fSimData.fSimulationSettings.filterBoundaryEqs = true;
    sim.fSimData.fSimulationSettings.axiSymmetricSimulation = false;
    sim.fSimData.fSimulationSettings.nThreads = 8;


    sim.fSimData.fOutputSettings.resultsDir = "results/concentricTime/";
    sim.fSimData.fOutputSettings.printGmeshVtk = true;
    sim.fSimData.fOutputSettings.printGmeshTxt = true;
    sim.fSimData.fOutputSettings.printCmeshTxt = true;
    sim.fSimData.fOutputSettings.printCmeshVtk = true;
    sim.fSimData.fOutputSettings.vtkSol = true;
    sim.fSimData.fOutputSettings.vtkResolution = 0;

}

void ConfigureFreqConcentricCase(TPZAcousticsSimulation &sim){
    ConfigureConcentricCase(sim);
    sim.fSimData.fSimulationSettings.simType = SPZAcousticData::ESimulationType::frequencyDomain;
    sim.fSimData.fOutputSettings.resultsDir = "results/concentricFreq/";
    sim.fSimData.fFourierSettings.nSamples = 200;
    sim.fSimData.fFourierSettings.wMax = sim.fSimData.fSourceSettings.centralFrequency * 3;
    sim.fSimData.fFourierSettings.alphaFreqShift = -1;
    sim.fSimData.fSimulationSettings.nTimeSteps = 361;
    sim.fSimData.fSimulationSettings.isCflBound = false;
}

void ConfigureHomogeneousCase(TPZAcousticsSimulation &sim){
    sim.fSimData.fSourceSettings.posX = 0;
    sim.fSimData.fSourceSettings.posY = 0;
    sim.fSimData.fSourceSettings.amplitude = 1;
    sim.fSimData.fSourceSettings.peakTime = 0.01;
    sim.fSimData.fSourceSettings.centralFrequency = 2*M_PI*100;

    sim.fSimData.fSimulationSettings.simType = SPZAcousticData::ESimulationType::timeDomain;
    sim.fSimData.fSimulationSettings.meshName = "wellMesh.geo";
    sim.fSimData.fSimulationSettings.boundType.Resize(1);
    sim.fSimData.fSimulationSettings.boundType[0] = SPZAcousticData::EBoundType::softwall;
    sim.fSimData.fSimulationSettings.nElemPerLambda = 10;
    sim.fSimData.fSimulationSettings.totalTime = 12 * sim.fSimData.fSourceSettings.peakTime;
    sim.fSimData.fSimulationSettings.cfl = 0.2;
    sim.fSimData.fSimulationSettings.nTimeSteps = 100;
    sim.fSimData.fSimulationSettings.isCflBound = true;
    sim.fSimData.fSimulationSettings.pOrder = 2;
    sim.fSimData.fSimulationSettings.filterBoundaryEqs = true;
    sim.fSimData.fSimulationSettings.axiSymmetricSimulation = false;
    sim.fSimData.fSimulationSettings.nThreads = 8;


    sim.fSimData.fOutputSettings.resultsDir = "results/homogeneousTime/";
    sim.fSimData.fOutputSettings.printGmeshVtk = true;
    sim.fSimData.fOutputSettings.printGmeshTxt = true;
    sim.fSimData.fOutputSettings.printCmeshTxt = true;
    sim.fSimData.fOutputSettings.printCmeshVtk = true;
    sim.fSimData.fOutputSettings.vtkSol = true;
    sim.fSimData.fOutputSettings.vtkResolution = 0;

}

void ConfigureFreqHomogeneousCase(TPZAcousticsSimulation &sim){
    ConfigureHomogeneousCase(sim);

    sim.fSimData.fSimulationSettings.nTimeSteps = 300;
    sim.fSimData.fSimulationSettings.isCflBound = true;

    sim.fSimData.fSimulationSettings.simType = SPZAcousticData::ESimulationType::frequencyDomain;
    sim.fSimData.fOutputSettings.resultsDir = "results/homogeneousFreq/";
    sim.fSimData.fFourierSettings.nSamples = 150;
    sim.fSimData.fFourierSettings.wMax = sim.fSimData.fSourceSettings.centralFrequency * 3;
    sim.fSimData.fFourierSettings.alphaFreqShift = -1;
}

void ConfigureAxiHomogeneousCase(TPZAcousticsSimulation &sim){
    sim.fSimData.fSourceSettings.posX = 0;
    sim.fSimData.fSourceSettings.posY = 10;
    sim.fSimData.fSourceSettings.amplitude = 1;
    sim.fSimData.fSourceSettings.peakTime = 0.005;
    sim.fSimData.fSourceSettings.centralFrequency = 2*M_PI*200;

    sim.fSimData.fSimulationSettings.simType = SPZAcousticData::ESimulationType::timeDomain;
    sim.fSimData.fSimulationSettings.meshName = "axiSymmetricMesh.geo";
    sim.fSimData.fSimulationSettings.boundType.Resize(2);
    sim.fSimData.fSimulationSettings.boundType[0] = SPZAcousticData::EBoundType::hardwall;
    sim.fSimData.fSimulationSettings.boundType[1] = SPZAcousticData::EBoundType::softwall;
    sim.fSimData.fSimulationSettings.nElemPerLambda = 10;
    sim.fSimData.fSimulationSettings.totalTime = 12 * sim.fSimData.fSourceSettings.peakTime;
    sim.fSimData.fSimulationSettings.cfl = 0.2;
    sim.fSimData.fSimulationSettings.nTimeSteps = 100;
    sim.fSimData.fSimulationSettings.isCflBound = true;
    sim.fSimData.fSimulationSettings.pOrder = 2;
    sim.fSimData.fSimulationSettings.filterBoundaryEqs = true;
    sim.fSimData.fSimulationSettings.axiSymmetricSimulation = false;
    sim.fSimData.fSimulationSettings.nThreads = 8;


    sim.fSimData.fOutputSettings.resultsDir = "results/axiHomogeneousTime/";
    sim.fSimData.fOutputSettings.printGmeshVtk = true;
    sim.fSimData.fOutputSettings.printGmeshTxt = true;
    sim.fSimData.fOutputSettings.printCmeshTxt = true;
    sim.fSimData.fOutputSettings.printCmeshVtk = true;
    sim.fSimData.fOutputSettings.vtkSol = true;
    sim.fSimData.fOutputSettings.vtkResolution = 0;

}

void ConfigureFreqAxiHomogeneousCase(TPZAcousticsSimulation &sim){
    ConfigureAxiHomogeneousCase(sim);

    sim.fSimData.fSimulationSettings.nTimeSteps = 300;
    sim.fSimData.fSimulationSettings.isCflBound = true;

    sim.fSimData.fSimulationSettings.simType = SPZAcousticData::ESimulationType::frequencyDomain;
    sim.fSimData.fOutputSettings.resultsDir = "results/axiHomogeneousFreq/";
    sim.fSimData.fFourierSettings.nSamples = 150;
    sim.fSimData.fFourierSettings.wMax = sim.fSimData.fSourceSettings.centralFrequency * 3;
    sim.fSimData.fFourierSettings.alphaFreqShift = -1;
}