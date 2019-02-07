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
    TPZAcousticAnalysis(TPZAcousticCompMesher * compMesher, const int &nThreads);
    virtual void InitializeComputations() = 0;
    virtual void SetUpGaussianSource(const REAL & wZero, const REAL & peakTime, const REAL & amplitude) = 0;
    virtual void RunSimulationSteps(const REAL &totalTime, const int &nTimeSteps) = 0;

    void PostProcess(int vtkResolution, std::string &prefix);
protected:


    virtual void InitializeSolver() = 0;

    const bool fFilterBoundaryEquations;
    const int fNThreads;
    int fNTimeSteps;
    TPZAnalysis fPzAnalysis;
    TPZStepSolver<STATE> fStepSolver;

    TPZAcousticCompMesher *fCompMesher;

    TPZFMatrix<STATE> fProbeSolution;
    TPZFMatrix<STATE> fTimeDomainSolution;
};


#endif