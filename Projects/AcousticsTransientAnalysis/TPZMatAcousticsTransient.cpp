//
// Created by Francisco Teixeira Orlandini on 3/6/18.
//

#include <pzaxestools.h>
#include <pzvec_extras.h>
#include "TPZMatAcousticsTransient.h"
#include "pzbndcond.h"

TPZMatAcousticsTransient::TPZMatAcousticsTransient() : TPZDiscontinuousGalerkin(), fRho(-1), fVelocity(-1){

}
TPZMatAcousticsTransient::TPZMatAcousticsTransient(int id) : TPZDiscontinuousGalerkin(id), fRho(-1), fVelocity(-1){

}
TPZMatAcousticsTransient::TPZMatAcousticsTransient(int id, const REAL &rho, const REAL &velocity) :
        TPZDiscontinuousGalerkin(id), fRho(rho), fVelocity(velocity){

}
TPZMatAcousticsTransient::TPZMatAcousticsTransient(const TPZMatAcousticsTransient &mat) : TPZDiscontinuousGalerkin(mat),
                                                                     fRho(mat.fRho), fVelocity(mat.fVelocity){

}
TPZMatAcousticsTransient::~TPZMatAcousticsTransient(){

}

void TPZMatAcousticsTransient::FillDataRequirements(TPZMaterialData &data)
{
    data.SetAllRequirements(false);
    data.fNeedsSol = true;
}

void TPZMatAcousticsTransient::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
    TPZFMatrix<REAL> &phi = data.phi;
    if(data.HSize == 0) return;
    const int nshape = phi.Rows();
    for(int i = 0; i < nshape; i++){
        for(int j = 0; j < nshape; j++){
            ek(i, j) += weight*data.phi(i,0)*data.phi(j,0)/(fRho*fVelocity*fVelocity*fDeltaT*fDeltaT);
            ek(i, j) -= weight*data.dphi(0,i)*data.dphi(0,j)/fRho;
            ek(i, j) -= weight*data.dphi(1,i)*data.dphi(1,j)/fRho;
        }//for j
    }//for i
}

void TPZMatAcousticsTransient::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef) {
    TPZFMatrix<REAL> &phi = data.phi;
    const int nshape = phi.Rows();
    STATE sourceVal = 0;
    STATE prevSol = data.sol[0][0];
    STATE prevPrevSol = data.sol[1][0];
    fSource(fCurrentTime,sourceVal);
    for(int i = 0; i < nshape; i++){
        ef(i,0) += weight * phi(i,0) * sourceVal;
        ef(i,0) += weight * phi(i,0) * 2* prevSol/(fRho * fVelocity * fVelocity * fDeltaT * fDeltaT);
        ef(i,0) -= weight * phi(i,0) * 1* prevPrevSol/(fRho * fVelocity * fVelocity * fDeltaT * fDeltaT);
    }//for i

}

void TPZMatAcousticsTransient::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek,
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

void TPZMatAcousticsTransient::SetSource(const std::function<void(const REAL &time, STATE &val)> &source) {
    fSource = source;
}

void TPZMatAcousticsTransient::SetCurrentTime(const REAL &time) {
    fCurrentTime = time;
}

void TPZMatAcousticsTransient::SetDeltaT(const REAL &deltaT) {
    fDeltaT = deltaT;
}

int TPZMatAcousticsTransient::VariableIndex(const std::string &name) {
    if(!strcmp("Pressure",name.c_str()))    return  1;
    return -1;
}

int TPZMatAcousticsTransient::NSolutionVariables(int var) {
    switch(var){
        case 1://pressure
            return 1;
        default:
            DebugStop();
    }
    return 0;
}
void TPZMatAcousticsTransient::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout) {
    Solout.Resize( this->NSolutionVariables(var));


    switch(var){
        case 1://pressure
            Solout = data.sol[0];
            break;
        default:
            DebugStop();
    }
}

int TPZMatAcousticsTransient::Dimension() const {
    return 2;
}

int TPZMatAcousticsTransient::NStateVariables() {
    return 1;
}

void TPZMatAcousticsTransient::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataL, REAL weight,
                                              TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) {
    DebugStop();
}
