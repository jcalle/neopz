//
// Created by Francisco Teixeira Orlandini on 3/6/18.
//

#include <pzaxestools.h>
#include <pzvec_extras.h>
#include "TPZMatAcousticsPml.h"

TPZMatAcousticsPml::TPZMatAcousticsPml(const int id,const TPZMatAcousticsFourier &mat,
                                       const bool &att_x, const REAL &pmlBeginX,
                                       const bool &att_y, const REAL &pmlBeginY,
                                       const REAL &alphaMax, const REAL &d) :
        TPZMatAcousticsFourier(mat), fAttX (att_x), fAttY (att_y),
        fPmlBeginX (pmlBeginX), fPmlBeginY (pmlBeginY),
        fAlphaMax (alphaMax), fD (d) , fW(-1)
{
    this->SetId(id);
    if(fAlphaMax < 0) DebugStop(); //for the attenuation to happen
                                   // this value must be positive
    if(fD < 0) DebugStop(); // pml width must be positive
    if(!fAttX && !fAttY) DebugStop();//a pml that doesnt attenuate
                                     // in any direction?
}

TPZMatAcousticsPml::~TPZMatAcousticsPml(){
}

void TPZMatAcousticsPml::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    TPZVec<REAL> &x = data.x;
    STATE dx =1., dy =1., d0;
    d0 = 3. * fAlphaMax * fVelocity / (2 * fD);

    if(fAttX){
        dx = 1. - SPZAlwaysComplex<STATE>::type(0,1) * (1./fW) * d0 * ((x[0]-fPmlBeginX) / fD)*((x[0]-fPmlBeginX) / fD);
    }
    if(fAttY){
        dy = 1. - SPZAlwaysComplex<STATE>::type(0,1) * (1./fW) * d0 * ((x[1]-fPmlBeginY) / fD)*((x[1]-fPmlBeginY) / fD);
    }
    TPZFMatrix<REAL> &phi = data.phi;
    const int nshape = phi.Rows();

    for(int i = 0; i < nshape; i++){
        for(int j = 0; j < nshape; j++){
            ek(i, j) += -1. * fW * fW * weight*data.phi(i,0)*data.phi(j,0)/(fRho * fVelocity * fVelocity);

            ek(i, j) += weight*data.dphix(0,i)*data.dphix(0,j)/(fRho * dx * dx);
            ek(i, j) += weight*data.dphix(1,i)*data.dphix(1,j)/(fRho * dy * dy);
        }//for j
    }//for i
}

int TPZMatAcousticsPml::IntegrationRuleOrder(TPZVec<int> &elPMaxOrder) const
{
    int pmax = 0;
    for (int ip=0;  ip<elPMaxOrder.size(); ip++)
    {
        if(elPMaxOrder[ip] > pmax) pmax = elPMaxOrder[ip];
    }
    const int integrationorder = 4+2*pmax;

    return  integrationorder;
}

STATE TPZMatAcousticsPml::GetW() const {
    return fW;
}

void TPZMatAcousticsPml::SetW(STATE fW) {
    TPZMatAcousticsPml::fW = fW;
}
