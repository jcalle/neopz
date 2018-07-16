//
// Created by Francisco Teixeira Orlandini on 3/6/18.
//

#ifndef PZ_TPZMATACOUSTICSH_H
#define PZ_TPZMATACOUSTICSH_H

#include <pzdiscgal.h>

class TPZMatAcousticsH1 : public TPZDiscontinuousGalerkin {
public:
    TPZMatAcousticsH1();//this does not exist
    explicit TPZMatAcousticsH1(int id);
    TPZMatAcousticsH1(int id, const REAL &rho, const REAL &velocity);
    TPZMatAcousticsH1(const TPZMatAcousticsH1 &mat);
    ~TPZMatAcousticsH1() override;
    void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
    void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) override;
    void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataL, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) override;
    int Dimension() const override;
    int NStateVariables() override;
    int VariableIndex(const std::string &name) override;
    int NSolutionVariables(int var) override;
    void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout) override;
    void SetExactSol(void (*exactSol)(const TPZVec<REAL> &, TPZVec<STATE> &, TPZFMatrix<STATE> &));
private:
    REAL fRho;
    REAL fVelocity;
    void (*fExactSol)(const TPZVec<REAL> &coord, TPZVec<STATE> &result,
                  TPZFMatrix<STATE> &grad);
};


#endif //PZ_TPZMATACOUSTICSH_H
