//
// Created by Francisco Teixeira Orlandini on 3/6/18.
//

#ifndef PZ_TPZMATACOUSTICSHTRANSIENT_H
#define PZ_TPZMATACOUSTICSHTRANSIENT_H

#include <pzdiscgal.h>
#include <functional>

class TPZMatAcousticsTransient : public TPZDiscontinuousGalerkin {
public:
    TPZMatAcousticsTransient();//this does not exist
    explicit TPZMatAcousticsTransient(int id);
    TPZMatAcousticsTransient(int id, const REAL &rho, const REAL &velocity);
    TPZMatAcousticsTransient(const TPZMatAcousticsTransient &mat);
    ~TPZMatAcousticsTransient() override;
    void FillDataRequirements(TPZMaterialData &data) override;
    void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
    void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef) override;
    void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) override;
    void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataL, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) override;
    int Dimension() const override;
    int NStateVariables() override;
    int VariableIndex(const std::string &name) override;
    int NSolutionVariables(int var) override;
    void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout) override;
    void SetSource(const std::function<void (const REAL &time, STATE &val)> &source);
    void SetCurrentTime(const REAL& time);
    void SetDeltaT(const REAL& deltat);
protected:
    REAL fRho;
    REAL fVelocity;
    REAL fDeltaT;
    REAL fCurrentTime;
    REAL fNewmarkBeta =0.25;
    std::function<void (const REAL &time, STATE &val)> fSource;
};


#endif //PZ_TPZMATACOUSTICSHTRANSIENT_H
