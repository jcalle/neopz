#ifndef TPZACOUSTICFREQDOMAINANALYSIS_H
#define TPZACOUSTICFREQDOMAINANALYSIS_H

#include "TPZAcousticAnalysis.h"

class TPZAcousticFreqDomainAnalysis : public TPZAcousticAnalysis{
public:

    /**
     * This constructor will not be generated
     */
    TPZAcousticFreqDomainAnalysis() = delete;
    /**
     * This is the only constructor that will be generated.
     * @param compMesher The instance of TPZAcousticCompMesher that will be used for the analysis
     */
    TPZAcousticFreqDomainAnalysis(TPZAcousticCompMesher * compMesher, const int &nThreads);
    void InitializeComputations() final;
    void SetUpGaussianSource(const REAL & wZero, const REAL & peakTime, const REAL & amplitude) final;


    void SetUpFourierSettings(const REAL &wMax, const int &nSamples, const REAL &aFactorUser);
    void RunSimulationSteps(const REAL &totalTime, const int &nTimeSteps) final;

protected:
    std::function<void (const STATE &currentW, STATE &val)> fFreqDomainSource;
    REAL fAFactor;
    REAL fWMax;
    int fNSamples;
    bool fAmIready;
};

#endif