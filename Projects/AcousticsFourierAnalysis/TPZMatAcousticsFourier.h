//
// Created by Francisco Teixeira Orlandini on 3/6/18.
//

#ifndef PZ_TPZMATACOUSTICSH_H
#define PZ_TPZMATACOUSTICSH_H

#include <pzdiscgal.h>
#include <functional>

class TPZMatAcousticsFourier : public TPZDiscontinuousGalerkin {
public:
    enum EWhichMatrix {M = 0, K, NDefined};
    TPZMatAcousticsFourier();//this does not exist
    explicit TPZMatAcousticsFourier(int id);
    TPZMatAcousticsFourier(int id, const REAL &rho, const REAL &velocity);
    TPZMatAcousticsFourier(const TPZMatAcousticsFourier &mat);
    ~TPZMatAcousticsFourier() override;
    void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
    void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef) override;
    void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) override;
    void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataL, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) override;
    int Dimension() const override;
    int NStateVariables() override;
    int VariableIndex(const std::string &name) override;
    int NSolutionVariables(int var) override;
    void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout) override;
    void SetExactSol(void (*exactSol)(const TPZVec<REAL> &, TPZVec<STATE> &, TPZFMatrix<STATE> &));

    void SetSourceFunc(const STATE &fSourceFunc);

    void SetAssemblingMatrix(EWhichMatrix mat);
protected:
    REAL fRho;
    REAL fVelocity;
    void (*fExactSol)(const TPZVec<REAL> &coord, TPZVec<STATE> &result,
                  TPZFMatrix<STATE> &grad);
    STATE fSourceFunc;

private:
    EWhichMatrix fAssembling;
};


#endif //PZ_TPZMATACOUSTICSH_H
