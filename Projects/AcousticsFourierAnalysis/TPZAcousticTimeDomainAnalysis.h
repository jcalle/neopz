#ifndef TPZACOUSTICTIMEDOMAINANALYSIS_H
#define TPZACOUSTICTIMEDOMAINANALYSIS_H

#include "TPZAcousticAnalysis.h"

class TPZAcousticTimeDomainAnalysis : public TPZAcousticAnalysis{
public:

    /**
     * This constructor will not be generated
     */
    TPZAcousticTimeDomainAnalysis() = delete;
    /**
     * This is the only constructor that will be generated.
     * @param compMesher The instance of TPZAcousticCompMesher that will be used for the analysis
     */
    TPZAcousticTimeDomainAnalysis(TPZAcousticCompMesher * compMesher, const int &nThreads, const bool &filter = true);
    void InitializeComputations() final;
    void SetUpGaussianSource(const REAL & wZero, const REAL & peakTime, const REAL & amplitude) final;

    void RunSimulationSteps(const REAL &totalTime, const int &nTimeSteps) final;

protected:

    void InitializeSolver() final;
    std::function<void (const REAL &time, STATE &val)> fTimeDomainSource;
    bool fAmIready;
};

#endif