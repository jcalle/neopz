//
// Created by Francisco Teixeira Orlandini on 3/6/18.
//

#include <pzaxestools.h>
#include <pzvec_extras.h>
#include "TPZMatAcousticsH1.h"

TPZMatAcousticsH1::TPZMatAcousticsH1() : TPZDiscontinuousGalerkin(), fRho(-1), fVelocity(-1){

}
TPZMatAcousticsH1::TPZMatAcousticsH1(int id) : TPZDiscontinuousGalerkin(id), fRho(-1), fVelocity(-1){

}
TPZMatAcousticsH1::TPZMatAcousticsH1(int id, const REAL &rho, const REAL &velocity) :
        TPZDiscontinuousGalerkin(id), fRho(rho), fVelocity(velocity){

}
TPZMatAcousticsH1::TPZMatAcousticsH1(const TPZMatAcousticsH1 &mat) : TPZDiscontinuousGalerkin(mat),
                                                                     fRho(mat.fRho), fVelocity(mat.fVelocity){

}
TPZMatAcousticsH1::~TPZMatAcousticsH1(){

}
void TPZMatAcousticsH1::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {

}

void TPZMatAcousticsH1::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek,
                                     TPZFMatrix<STATE> &ef, TPZBndCond &bc){

}
void TPZMatAcousticsH1::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout) {

}
int TPZMatAcousticsH1::IntegrationRuleOrder(TPZVec<int> &elPMaxOrder) const {

}