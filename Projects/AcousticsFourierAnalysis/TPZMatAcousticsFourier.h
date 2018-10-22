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

    void SetSource(const std::function<void (const STATE &w, STATE &val)> &source);

    void SetAssemblingMatrix(EWhichMatrix mat);
    STATE GetW() const;

    void SetW(STATE fW);

protected:
    STATE fW;
    REAL fRho;
    REAL fVelocity;
    std::function<void (const STATE &time, STATE &val)> fSource;

private:
    EWhichMatrix fAssembling;
};


#endif //PZ_TPZMATACOUSTICSH_H
