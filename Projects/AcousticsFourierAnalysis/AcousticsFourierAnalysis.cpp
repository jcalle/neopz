#include "TPZAcousticsSimulation.h"

enum ESimulationCases{ freqAxiSymmetric, freqHomogeneous2D, freqHeterogeneous2D,
                      timeAxiSymmetric, timeHomogeneous2D, timeHeterogeneous2D};

void ConfigureHomogeneousCase(TPZAcousticsSimulation &sim);
void ConfigureFreqHomogeneousCase(TPZAcousticsSimulation &sim);

int main(int argc, char *argv[]) {
    TPZAcousticsSimulation sim;
    ESimulationCases simCase = ESimulationCases::freqHomogeneous2D;
    switch(simCase){
        case ESimulationCases::freqAxiSymmetric:
            break;
        case ESimulationCases::freqHomogeneous2D:
            ConfigureFreqHomogeneousCase(sim);
            break;
        case ESimulationCases::freqHeterogeneous2D:
            break;
        case ESimulationCases::timeAxiSymmetric:
            break;
        case ESimulationCases::timeHomogeneous2D:
            break;
        case ESimulationCases::timeHeterogeneous2D:
            break;
        default:
            PZError<<"What are you trying to run?"<<std::endl;
            DebugStop();
    }
    sim.RunSimulation();
}


void ConfigureHomogeneousCase(TPZAcousticsSimulation &sim){

    sim.fSimData.fSimulationSettings.simType = SPZAcousticData::ESimulationType::timeDomain;
    sim.fSimData.fSimulationSettings.meshName = "wellMesh.geo";
    sim.fSimData.fSimulationSettings.boundType = SPZAcousticData::EBoundType::softwall;
    sim.fSimData.fSimulationSettings.nElemPerLambda = 12;
    sim.fSimData.fSimulationSettings.totalTime = 0.1;
    sim.fSimData.fSimulationSettings.cfl = 0.2;
    sim.fSimData.fSimulationSettings.nTimeSteps = 100;
    sim.fSimData.fSimulationSettings.isCflBound = true;
    sim.fSimData.fSimulationSettings.pOrder = 2;
    sim.fSimData.fSimulationSettings.filterBoundaryEqs = true;
    sim.fSimData.fSimulationSettings.axiSymmetricSimulation = false;
    sim.fSimData.fSimulationSettings.nThreads = 8;


    sim.fSimData.fSourceSettings.posX = 0;
    sim.fSimData.fSourceSettings.posY = 0;
    sim.fSimData.fSourceSettings.amplitude = 1;
    sim.fSimData.fSourceSettings.peakTime = 0.01;
    sim.fSimData.fSourceSettings.centralFrequency = 100;


    sim.fSimData.fOutputSettings.resultsDir = "results/";
    sim.fSimData.fOutputSettings.printGmeshVtk = true;
    sim.fSimData.fOutputSettings.printGmeshTxt = true;
    sim.fSimData.fOutputSettings.printCmeshTxt = true;
    sim.fSimData.fOutputSettings.printCmeshVtk = true;
    sim.fSimData.fOutputSettings.vtkSol = true;
    sim.fSimData.fOutputSettings.vtkResolution = 0;

}

void ConfigureFreqHomogeneousCase(TPZAcousticsSimulation &sim){
    ConfigureHomogeneousCase(sim);
    sim.fSimData.fSimulationSettings.simType = SPZAcousticData::ESimulationType::frequencyDomain;
    sim.fSimData.fFourierSettings.nSamples = 50;
    sim.fSimData.fFourierSettings.wMax = 300;
    sim.fSimData.fFourierSettings.alphaFreqShift = -1;
}