//
// Created by Francisco Teixeira Orlandini on 3/6/18.
//

#include <pzaxestools.h>
#include <pzvec_extras.h>
#include "TPZMatAcousticsH1.h"
#include "pzbndcond.h"

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
    TPZFMatrix<REAL> &phi = data.phi;
    const int nshape = phi.Rows();
    TPZVec<STATE> sol;
    TPZFMatrix<STATE> grad;
    (*fExactSol)(data.x,sol,grad);
    for(int i = 0; i < nshape; i++){
        for(int j = 0; j < nshape; j++){
            ek(i, j) += weight*data.phi(i,0)*data.phi(j,0);
        }//for j
        ef(i,0) += (STATE)weight*(STATE)data.phi(i,0)*sol[0];
    }//for i
}

void TPZMatAcousticsH1::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek,
                                     TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    return;
    TPZFMatrix<REAL> &phi = data.phi;
    const int phr = phi.Rows();
    int in, jn;

    switch (bc.Type()){

        // Dirichlet condition
        case 0 : {
            for(in = 0 ; in < phr; in++) {
                ef(in,0) += (STATE)TPZMaterial::gBigNumber * bc.Val2()(0,0) * (STATE)phi(in,0) * (STATE)weight;
                for (jn = 0 ; jn < phr; jn++) {
                    ek(in,jn) += TPZMaterial::gBigNumber * phi(in,0) * phi(jn,0) * weight;
                }//jn
            }//in
            break;
        }

            // Neumann condition
        case 1 : {
            for(in = 0 ; in < phr; in++) {
                ef(in,0) += bc.Val2()(0,0) * (STATE)phi(in,0) * (STATE)weight;
            }//in
            break;
        }

        default:{
            std::cout << __PRETTY_FUNCTION__ << " at line " << __LINE__ << " not implemented\n";
        }
    }//switch

}

int TPZMatAcousticsH1::VariableIndex(const std::string &name) {
    if(!strcmp("Pressure",name.c_str()))    return  1;
}

int TPZMatAcousticsH1::NSolutionVariables(int var) {
    switch(var){
        case 1://pressure
            return 1;
        default:
            DebugStop();
    }
    return 0;
}
void TPZMatAcousticsH1::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout) {
    Solout.Resize( this->NSolutionVariables(var));


    switch(var){
        case 1://pressure
            Solout = data.sol[0];
            for(int i = 0; i< Solout.size(); i++){
                Solout[i] = std::real(Solout[i]);
            }
            break;
        default:
            DebugStop();
    }
}

int TPZMatAcousticsH1::Dimension() const {
    return 2;
}

int TPZMatAcousticsH1::NStateVariables() {
    return 1;
}

void TPZMatAcousticsH1::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataL, REAL weight,
                                              TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) {
    DebugStop();
}

void TPZMatAcousticsH1::SetExactSol(void (*exactSol)(const TPZVec<REAL> &, TPZVec<STATE> &, TPZFMatrix<STATE> &)) {
    TPZMatAcousticsH1::fExactSol = exactSol;
}
