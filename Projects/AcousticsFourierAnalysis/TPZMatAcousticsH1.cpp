//
// Created by Francisco Teixeira Orlandini on 3/6/18.
//

#include <pzaxestools.h>
#include <pzvec_extras.h>
#include "TPZMatAcousticsH1.h"
#include "pzbndcond.h"

TPZMatAcousticsH1::TPZMatAcousticsH1() : TPZDiscontinuousGalerkin(), fRho(-1), fVelocity(-1), fAssembling(NDefined){

}
TPZMatAcousticsH1::TPZMatAcousticsH1(int id) : TPZDiscontinuousGalerkin(id), fRho(-1), fVelocity(-1),
                                               fAssembling(NDefined){

}
TPZMatAcousticsH1::TPZMatAcousticsH1(int id, const REAL &rho, const REAL &velocity) :
        TPZDiscontinuousGalerkin(id), fRho(rho), fVelocity(velocity) , fAssembling(NDefined){

}
TPZMatAcousticsH1::TPZMatAcousticsH1(const TPZMatAcousticsH1 &mat) : TPZDiscontinuousGalerkin(mat),
                                                                     fRho(mat.fRho), fVelocity(mat.fVelocity),
                                                                     fAssembling(NDefined){

}
TPZMatAcousticsH1::~TPZMatAcousticsH1(){

}

void TPZMatAcousticsH1::SetAssemblingMatrix(TPZMatAcousticsH1::EWhichMatrix mat) {
    fAssembling = mat;
}
void TPZMatAcousticsH1::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
    TPZFMatrix<REAL> &phi = data.phi;
    const int nshape = phi.Rows();

    switch(fAssembling){
        case M:
            for(int i = 0; i < nshape; i++){
                for(int j = 0; j < nshape; j++){
                    ek(i, j) += weight*data.phi(i,0)*data.phi(j,0)/(fRho * fVelocity * fVelocity);
                }//for j
            }//for i
            break;
        case K:
            for(int i = 0; i < nshape; i++){
                for(int j = 0; j < nshape; j++){
                    ek(i, j) += weight*data.dphi(0,i)*data.dphi(0,j)/fRho;
                    ek(i, j) += weight*data.dphi(1,i)*data.dphi(1,j)/fRho;
                }//for j
            }//for i
            break;
        case NDefined:
            DebugStop();
    }
}

void TPZMatAcousticsH1::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef) {
    TPZFMatrix<REAL> &phi = data.phi;
    const int nshape = phi.Rows();

    for(int i = 0; i < nshape; i++){
        ef(i,0) += weight * phi(i,0) * fSourceFunc / (fRho* fVelocity * fVelocity);
    }//for i

}

void TPZMatAcousticsH1::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek,
                                     TPZFMatrix<STATE> &ef, TPZBndCond &bc){
//    return;//Neumann 0 you do nothing, Dirichlet 0 will be filtered
    TPZFMatrix<REAL> &phi = data.phi;
    const int phr = phi.Rows();
    int in, jn;

    switch (bc.Type()){

        // Dirichlet condition
        case 0 : {
            for(in = 0 ; in < phr; in++) {
                ef(in,0) += (STATE)TPZMaterial::gBigNumber * bc.Val2()(0,0) * (STATE)phi(in,0) * (STATE)weight;
                for (jn = 0 ; jn < phr; jn++) {
                    ek(in,jn) += weight * TPZMaterial::gBigNumber * phi(in,0) * phi(jn,0)  * (1/fRho);
                }//jn
            }//in
            break;
        }

            // Neumann condition
        case 1 : {
            for(in = 0 ; in < phr; in++) {
                ef(in,0) +=(STATE)weight * bc.Val2()(0,0) * (STATE)phi(in,0) * (1/fRho);
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
    return -1;
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

void TPZMatAcousticsH1::SetSourceFunc(const STATE &sourceFunc) {
    fSourceFunc = sourceFunc;
}
