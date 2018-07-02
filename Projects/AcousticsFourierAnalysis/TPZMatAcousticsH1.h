//
// Created by Francisco Teixeira Orlandini on 3/6/18.
//

#ifndef PZ_TPZMATACOUSTICSH!_H
#define PZ_TPZMATACOUSTICSH!_H

#include <pzdiscgal.h>

class TPZMatAcousticsH1 : public TPZDiscontinuousGalerkin {
public:
    TPZMatAcousticsH1();//this does not exist
    TPZMatAcousticsH1(int id);
    TPZMatAcousticsH1(int id, const REAL &rho, const REAL &velocity);
    TPZMatAcousticsH1(const TPZMatAcousticsH1 &mat);
    ~TPZMatAcousticsH1();
    void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
    void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) override;
    void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout) override;
    int IntegrationRuleOrder(TPZVec<int> &elPMaxOrder) const override;
private:
    REAL fRho;
    REAL fVelocity;
};


#endif //PZ_TPZMATACOUSTICSH!_H
