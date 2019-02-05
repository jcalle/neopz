#include <pzgmesh.h>
#include <pzvec.h>

struct SPZAcousticData{
//    enum pmltype{ TODO: tarefa pra 2019 (2020?)
//        xp=0,yp,xm,ym,xpyp,xmyp,xmym,xpym
//    };
    enum EBoundType{
        softwall = 0, hardwall = 1
    };
    enum ESimulationType{
        timeDomain = 0, frequencyDomain = 1
    };
    struct SPZSimulationSettings{
    	ESimulationType  simType;
        std::string meshName = "iDontExist.geo";
        TPZVec<SPZAcousticData::EBoundType> boundType;
        int nElemPerLambda;
        REAL totalTime;
        REAL cfl;
        int nTimeSteps;
        bool isCflBound;
        int pOrder;
        bool filterBoundaryEqs;
        bool axiSymmetricSimulation;
        int nThreads;
    };
    struct SPZSourceSettings{
        REAL posX;
        REAL posY;
        REAL amplitude;
        REAL peakTime;
        REAL centralFrequency;
    };
    struct SPZOutputSettings{
		std::string resultsDir = "results/";
        bool printGmeshVtk;
        bool printGmeshTxt;
        bool printCmeshTxt;
        bool printCmeshVtk;
        bool vtkSol;
        int vtkResolution;
    };

	struct SPZAcousticFourierData{
		int nSamples;
		REAL wMax;
		REAL alphaFreqShift;
	};

    SPZSimulationSettings fSimulationSettings;
    SPZSourceSettings fSourceSettings;
    SPZOutputSettings fOutputSettings;
    SPZAcousticFourierData fFourierSettings;
};

//struct SPZAcousticTransientData{
//
//};

class TPZAcousticsSimulation {
public:
	void RunSimulation();
	SPZAcousticData fSimData;
protected:
	REAL CalculateDeltaT(std::map<int,REAL> elSizeMap, std::map<int,REAL> velocityMap);
};
//
//class TPZAcousticsFourierSimulation : public TPZAcousticsSimulation{
//protected:
//	SPZAcousticFourierData fFourierData;
//
//};
//
//class TPZAcousticsTransientSimulation : public TPZAcousticsSimulation{
//	SPZAcousticTransientData fTransientData;
//
//};