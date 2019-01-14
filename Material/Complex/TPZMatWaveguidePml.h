//
// Created by Francisco Teixeira Orlandini on 3/6/18.
//

#ifndef PZ_TPZMATWAVEGUIDEPML_H
#define PZ_TPZMATWAVEGUIDEPML_H


#include "TPZMatModalAnalysis.h"

class TPZMatWaveguidePml : public virtual TPZMatModalAnalysis {
public:
    TPZMatWaveguidePml(const int id,const TPZMatModalAnalysis &mat,
                       const bool &att_x, REAL &pmlBeginX,
                       const bool &att_y, REAL &pmlBeginY,
                       const REAL &alphaMax, const REAL &d);
    TPZMatWaveguidePml(const TPZMatModalAnalysis &mat) = delete;//this does not exist
    TPZMatWaveguidePml(int id) = delete;//this does not exist
    TPZMatWaveguidePml() = delete;//this does not exist
    ~TPZMatWaveguidePml() = default;
    void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
    void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout) override;
    int IntegrationRuleOrder(TPZVec<int> &elPMaxOrder) const override;
protected:
    const bool fAttX;
    const REAL fPmlBeginX;
    const bool fAttY;
    const REAL fPmlBeginY;
    const REAL fAlphaMax;
    const REAL fD;
private:
    const REAL fTolAlpha = 1e-12;
};


#endif //PZ_TPZMATWAVEGUIDEPML_H
