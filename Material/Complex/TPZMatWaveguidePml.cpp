//
// Created by Francisco Teixeira Orlandini on 3/6/18.
//

#include <pzaxestools.h>
#include <pzvec_extras.h>
#include "TPZMatWaveguidePml.h"

TPZMatWaveguidePml::TPZMatWaveguidePml(const int id,const TPZMatModalAnalysis &mat,
                                       const bool &att_x, REAL &pmlBeginX,
                                       const bool &att_y, REAL &pmlBeginY,
                                       const REAL &alphaMax, const REAL &d) :
        TPZMatModalAnalysis(mat), fAttX (att_x), fAttY (att_y),
        fPmlBeginX (pmlBeginX), fPmlBeginY (pmlBeginY),
        fAlphaMax (alphaMax), fD (d)
{
    this->SetId(id);
    if(fAlphaMax < 0) DebugStop(); //for the attenuation to happen
                                   // this value must be positive
    if(fD < 0) DebugStop(); // pml width must be positive
    if(!fAttX && !fAttY) DebugStop();//a pml that doesnt attenuate
                                     // in any direction?
}

void TPZMatWaveguidePml::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    if(fAlphaMax < fTolAlpha){
        this->RealContribute(datavec,weight,ek,ef,fEr,fEr,fEr,fUr,fUr,fUr);
        return;
    }
    /*****************CALCULATE S PML PARAMETERS*************************
     * In the current application, the waveguide's cross section is always
     * in the xy-plane. Therefore, sz will always be unity, and omitted for
     * the folllowing calculations. The same principle applies, for instance,
     * for the z-component of the hcurl functions, the x and y components of
     * their curl and so on.
     */
    TPZManVector<REAL,3> x = datavec[ h1meshindex ].x;
    STATE sx = 1, sy = 1;
    if(fAttX){
        sx = 1. - imaginary * fAlphaMax * ((x[0]-fPmlBeginX) / fD ) * ((x[0]-fPmlBeginX) / fD );
    }
    if(fAttY){
        sy = 1. - imaginary * fAlphaMax * ((x[1]-fPmlBeginY) / fD ) * ((x[1]-fPmlBeginY) / fD );
    }
    const STATE uxx = fUr * sy / sx;
    const STATE uyy = fUr * sx / sy;
    const STATE uzz = fUr * sy * sx;
    const STATE exx = fEr * sy / sx;
    const STATE eyy = fEr * sx / sy;
    const STATE ezz = fEr * sy * sx;

    this->RealContribute(datavec,weight,ek,ef,exx,eyy,ezz,uxx,uyy,uzz);
}



void TPZMatWaveguidePml::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout)
{
    if(fAlphaMax < fTolAlpha){
        this->RealSolution(datavec,var,Solout,fEr,fEr,fUr,fUr);
        return;
    }
    /*****************CALCULATE S PML PARAMETERS*************************/
    TPZManVector<REAL,3> x = datavec[ h1meshindex ].x;
    REAL sx = 1, sy = 1;
    Solout.Resize(2);
    if(fAttX){
        sx = 1 + fAlphaMax * ((x[0]-fPmlBeginX) / fD ) * ((x[0]-fPmlBeginX) / fD );
    }
    if(fAttY){
        sy = 1 + fAlphaMax * ((x[1]-fPmlBeginY) / fD ) * ((x[1]-fPmlBeginY) / fD );
    }
    /**********************
     * The following is just a way to see the PML attenuation
     */
    const STATE uxx = fUr / sx;
    const STATE uyy = fUr / sy;
    const STATE exx = fEr / sx;
    const STATE eyy = fEr / sy;

    this->RealSolution(datavec,var,Solout,exx,eyy,uxx,uyy);
}

int TPZMatWaveguidePml::IntegrationRuleOrder(TPZVec<int> &elPMaxOrder) const
{
    if(fAlphaMax < fTolAlpha){
        return TPZMatModalAnalysis::IntegrationRuleOrder(elPMaxOrder);
    }
    int pmax = 0;
    for (int ip=0;  ip<elPMaxOrder.size(); ip++)
    {
        if(elPMaxOrder[ip] > pmax) pmax = elPMaxOrder[ip];
    }

    const int integrationorder = 4+2*pmax;

    return  integrationorder;
}
