//
// Created by Francisco Teixeira Orlandini on 3/6/18.
//

#include "TPZMatAcousticsPml.h"
#include <pzaxestools.h>
#include <pzvec_extras.h>
#include "TPZMatAcousticsFourier.h"
#include "pzbndcond.h"

TPZMatAcousticsFourier::TPZMatAcousticsFourier() : TPZDiscontinuousGalerkin(), fRho(-1), fVelocity(-1), fW(-1), fAssembling(NDefined){

}
TPZMatAcousticsFourier::TPZMatAcousticsFourier(int id) : TPZDiscontinuousGalerkin(id), fRho(-1), fVelocity(-1), fW(-1),
                                               fAssembling(NDefined){

}
TPZMatAcousticsFourier::TPZMatAcousticsFourier(int id, const REAL &rho, const REAL &velocity, const bool & isAxisymmetric) :
        TPZDiscontinuousGalerkin(id), fRho(rho), fVelocity(velocity) , fW(-1), fAssembling(NDefined),
        fIsAxisymmetric(isAxisymmetric){

}
TPZMatAcousticsFourier::TPZMatAcousticsFourier(const TPZMatAcousticsFourier &mat) : TPZDiscontinuousGalerkin(mat),
                                                                     fRho(mat.fRho),  fW(mat.fW),fVelocity(mat.fVelocity),
                                                                     fAssembling(NDefined){

}
TPZMatAcousticsFourier::~TPZMatAcousticsFourier(){

}

void TPZMatAcousticsFourier::SetAssemblingMatrix(TPZMatAcousticsFourier::EWhichMatrix mat) {
    fAssembling = mat;
}
void TPZMatAcousticsFourier::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
    REAL actualWeight = weight;
    if(fIsAxisymmetric){
        actualWeight *= data.x[0]*2*M_PI;
    }
    TPZFMatrix<REAL> &phi = data.phi;
    const int nshape = phi.Rows();

    switch(fAssembling){
        case M:
            for(int i = 0; i < nshape; i++){
                for(int j = 0; j < nshape; j++){
                    ek(i, j) += actualWeight*data.phi(i,0)*data.phi(j,0)/(fRho * fVelocity * fVelocity);
                }//for j
            }//for i
            break;
        case K:
            for(int i = 0; i < nshape; i++){
                for(int j = 0; j < nshape; j++){
                    ek(i, j) += actualWeight*data.dphix(0,i)*data.dphix(0,j)/fRho;
                    ek(i, j) += actualWeight*data.dphix(1,i)*data.dphix(1,j)/fRho;
                }//for j
            }//for i
            break;
        case NDefined:
            DebugStop();
    }
}

void TPZMatAcousticsFourier::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef) {
    TPZFMatrix<REAL> &phi = data.phi;
    const int nshape = phi.Rows();
    if(nshape > 1){
        DebugStop();
    }
    STATE sourceVal;
    fSource(GetW(),sourceVal);
    for(int i = 0; i < nshape; i++){
        ef(i,0) += phi(i,0) * sourceVal;
    }//for i

}

void TPZMatAcousticsFourier::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek,
                                     TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    REAL actualWeight = weight;
    if(fIsAxisymmetric){
        actualWeight *= data.x[0]*2*M_PI;
    }
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
                    ek(in,jn) += actualWeight * TPZMaterial::gBigNumber * phi(in,0) * phi(jn,0);
                }//jn
            }//in
            break;
        }

            // Neumann condition
        case 1 : {
            for(in = 0 ; in < phr; in++) {
                ef(in,0) +=(STATE)weight * bc.Val2()(0,0) * (STATE)phi(in,0);
            }//in
            break;
        }

        default:{
            std::cout << __PRETTY_FUNCTION__ << " at line " << __LINE__ << " not implemented\n";
        }
    }//switch

}

int TPZMatAcousticsFourier::VariableIndex(const std::string &name) {
    if(!strcmp("Pressure",name.c_str()))    return  1;
    return -1;
}

int TPZMatAcousticsFourier::NSolutionVariables(int var) {
    switch(var){
        case 1://pressure
            return 1;
        default:
            DebugStop();
    }
    return 0;
}
void TPZMatAcousticsFourier::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout) {
    Solout.Resize( this->NSolutionVariables(var));


    switch(var){
        case 1://pressure
            Solout = data.sol[0];
            break;
        default:
            DebugStop();
    }
}

int TPZMatAcousticsFourier::Dimension() const {
    return 2;
}

int TPZMatAcousticsFourier::NStateVariables() {
    return 1;
}

void TPZMatAcousticsFourier::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataL, REAL weight,
                                              TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) {
    DebugStop();
}

void TPZMatAcousticsFourier::SetSource(const std::function<void(const STATE &omega, STATE &val)> &source) {
    fSource = source;
}

STATE TPZMatAcousticsFourier::GetW() const {
    return fW;
}

void TPZMatAcousticsFourier::SetW(STATE fW) {
    TPZMatAcousticsFourier::fW = fW;
}

int TPZMatAcousticsFourier::IntegrationRuleOrder(int elPMaxOrder) const
{
    int order = 0;

    int pmax = elPMaxOrder;
    int integrationorder = 2*pmax;
    if (pmax < order) {
        integrationorder = order+pmax;
    }

//    if(this->fIsAxisymmetric){
//        if(this->)
//        order = fForcingFunction->PolynomialOrder();
//    }
    return  integrationorder;
}