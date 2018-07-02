//
// Created by Francisco Teixeira Orlandini on 3/6/18.
//

#ifndef PZ_TPZMATACOUSTICSPML_H
#define PZ_TPZMATACOUSTICSPML_H


#include "TPZMatAcousticsH1.h"

class TPZMatAcousticsPml : public TPZMatAcousticsH1 {
public:
    TPZMatAcousticsPml(const int id,const TPZMatAcousticsH1 &mat,
                       const bool &att_x, const REAL &pmlBeginX,
                       const bool &att_y, const REAL &pmlBeginY,
                       const REAL &alphaMax, const REAL &d);
    ~TPZMatAcousticsPml();
    void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
    void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout) override;
    int IntegrationRuleOrder(TPZVec<int> &elPMaxOrder) const override;
    STATE GetW() const;
    void SetW(STATE fW);
private:
    TPZMatAcousticsPml(const TPZMatAcousticsH1 &mat);//this does not exist
    TPZMatAcousticsPml(int id);//this does not exist
    TPZMatAcousticsPml();//this does not exist
    const bool fAttX;
    const bool fAttY;
    const REAL fPmlBeginX;
    const REAL fPmlBeginY;
    const REAL fAlphaMax;
    const REAL fD;
    STATE fW;
};


#endif //PZ_TPZMATACOUSTICSPML_H
