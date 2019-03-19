#ifndef TPZACOUSTICANALYSIS_H
#define TPZACOUSTICANALYSIS_H

#include "pzvec.h"
#include "pzstrmatrix.h"
#include "pzanalysis.h"
#include "pzysmp.h"
#include "pzstepsolver.h"

class TPZAcousticCompMesher;

/**
 * A class that creates a computational mesh based on a previously generated geometric mesh.
 * The geometric mesh must have been generated using a TPZAcousticGeoMesher in order to assure that all
 * relevant data structures have been properly filled.
 *
 */
class TPZAcousticAnalysis{
public:
    virtual ~TPZAcousticAnalysis() = default;
    /**
     * This constructor will not be generated
     */
    TPZAcousticAnalysis() = delete;
    /**
     * This is the only constructor that will be generated.
     * @param compMesher The instance of TPZAcousticCompMesher that will be used for the analysis
     */
    TPZAcousticAnalysis(TPZAcousticCompMesher * compMesher, const int &nThreads, const REAL &deltaT,
            const int &nTimeSteps,const bool &filter);
    virtual void InitializeComputations() = 0;
    virtual void SetUpGaussianSource(const REAL & wZero, const REAL & peakTime, const REAL & amplitude) = 0;
    virtual void RunSimulationSteps() = 0;
    void ProbeSolution(std::string &prefix);
    void PostProcess(int vtkResolution, std::string &prefix);
protected:
    void FindProbe(TPZCompEl *&probeComPEl, TPZGeoMesh *gmesh, const int &matIdSource);
    void FilterBoundaryEquations(TPZVec<int64_t> &activeEquations, int64_t &neq, int64_t &neqOriginal);
    virtual void InitializeSolver() = 0;

    bool isThereAProbe;

    TPZVec<int64_t > fActiveEquations;

    int64_t fNeqReduced;

    int64_t fNeqOriginal;

    bool fFilterBoundaryEquations;

    const int fNThreads;
    int fNTimeSteps;
    REAL fDeltaT;

    REAL fTotalTime;
    TPZAnalysis fPzAnalysis;
    TPZStepSolver<STATE> fStepSolver;

    TPZAcousticCompMesher *fCompMesher;

    TPZFMatrix<SPZAlwaysReal<STATE>::type> fProbeSolution;
    TPZFMatrix<STATE> fTimeDomainSolution;
};


#endif