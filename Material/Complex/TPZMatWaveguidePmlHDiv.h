//
// Created by Francisco Teixeira Orlandini on 3/6/18.
//

#ifndef PZ_TPZMATWAVEGUIDEPMLHDIV_H
#define PZ_TPZMATWAVEGUIDEPMLHDIV_H


#include "TPZMatWaveguidePml.h"
#include "TPZMatModalAnalysisHDiv.h"

class TPZMatWaveguidePmlHDiv : public TPZMatWaveguidePml, public TPZMatModalAnalysisHDiv {
public:
    TPZMatWaveguidePmlHDiv(const int id,const TPZMatModalAnalysis &mat,
                       const bool &att_x, REAL &pmlBeginX,
                       const bool &att_y, REAL &pmlBeginY,
                       const REAL &alphaMax, const REAL &d);
    ~TPZMatWaveguidePmlHDiv() = default;
    void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
    void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout) override;
private:
    TPZMatWaveguidePmlHDiv(const TPZMatModalAnalysis &mat);//this does not exist
    TPZMatWaveguidePmlHDiv(int id);//this does not exist
    TPZMatWaveguidePmlHDiv();//this does not exist
};


#endif //PZ_TPZMATWAVEGUIDEPMLHDIV_H
