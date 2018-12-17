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
        transient = 0, frequencyDomain = 1
    };
    struct SPZSimulationSettings{
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

    SPZSimulationSettings fSimulationSettings;
    SPZSourceSettings fSourceSettings;
    SPZOutputSettings fOutputSettings;
};

struct SPZAcousticFourierData{
	REAL wSample;
    REAL wMax;
    REAL alphaFreqShift;
};

struct SPZAcousticTransientData{

};

class TPZAcousticsSimulation {
public:

	void RunSimulation();

	virtual void SolveSystem() = 0;
	virtual void CreateCMesh(TPZCompMesh *&cmesh, TPZGeoMesh *gmesh, int pOrder, const std::string &prefix,
			const bool &print, const TPZVec<int> &matIdVec, const std::map<int,REAL> &rhoMap,
			const std::map<int,REAL> &velocityMap,
//          const REAL &alphaPml, const TPZVec< SPZAcousticData::pmltype > &pmlTypeVec,
							 const TPZVec<SPZAcousticData::EBoundType> &boundTypeVec) = 0;

	REAL CalculateDeltaT(const TPZVec<REAL> &elSizeVec, std::map<int,REAL> velocityMap);

	static void GetMaterialProperties(std::map<std::string,int> &materialNames, std::map<int,REAL> &rhoMap,
                           std::map<int,REAL> &velocityMap);

protected:

	void FilterBoundaryEquations(TPZCompMesh *cmesh, TPZVec<int64_t> &activeEquations, int &neq, int &neqOriginal);

	SPZAcousticData fSimData;
		
};

class TPZAcousticsFourierSimulation : public TPZAcousticsSimulation{
protected:
	SPZAcousticFourierData fFourierData;
	
};

class TPZAcousticsTransientSimulation : public TPZAcousticsSimulation{
	SPZAcousticTransientData fTransientData;
	
};