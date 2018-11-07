//
// Created by Francisco Teixeira Orlandini on 3/6/18.
//

#include <pzaxestools.h>
#include <pzvec_extras.h>
#include "TPZMatWaveguidePmlHDiv.h"

TPZMatWaveguidePmlHDiv::TPZMatWaveguidePmlHDiv(const int id,const TPZMatModalAnalysis &mat,
                                       const bool &att_x, REAL &pmlBeginX,
                                       const bool &att_y, REAL &pmlBeginY,
                                       const REAL &alphaMax, const REAL &d) :
        TPZMatWaveguidePml(id,mat,att_x,pmlBeginX, att_y,pmlBeginY,alphaMax, d) {
}

void TPZMatWaveguidePmlHDiv::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
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
    /*********************CREATE H1 FUNCTIONS****************************/
    TPZFNMatrix<12,REAL> phiH1 = datavec[ h1meshindex ].phi;
    TPZFNMatrix<36,REAL> dphiH1daxes = datavec[ h1meshindex ].dphix;
    TPZFNMatrix<3,REAL> dphiH1;
    TPZAxesTools<REAL>::Axes2XYZ(dphiH1daxes, dphiH1, datavec[ h1meshindex ].axes);
    TPZFNMatrix<3,REAL> gradPhiH1(phiH1.Rows() , 3 , 0.);
    for ( int iFunc = 0 ; iFunc < phiH1.Rows(); iFunc++ ) {

        gradPhiH1 ( iFunc , 0 ) = dphiH1 ( 0 , iFunc );
        gradPhiH1 ( iFunc , 1 ) = dphiH1 ( 1 , iFunc );
    }

    /*********************CREATE HDIV FUNCTIONS****************************/
    TPZFNMatrix<12,REAL> phiScaHCurl = datavec[ hcurlmeshindex ].phi;
    TPZManVector<REAL,3> xParametric = datavec[ h1meshindex ].xParametric;
    
    int phrq = datavec[ hcurlmeshindex ].fVecShapeIndex.NElements();
    //  std::cout<<"x"<<std::endl<<x[0]<<" "<<x[1]<<" "<<x[2]<<std::endl;
    
    TPZFNMatrix< 36 , REAL > phiVecHDiv(phrq , 3 , 0.);
    for (int iq = 0 ; iq < phrq ; iq++) {
        int ivecind = datavec[ hcurlmeshindex ].fVecShapeIndex[iq].first;
        int ishapeind = datavec[ hcurlmeshindex ].fVecShapeIndex[iq].second;
        
        phiVecHDiv(iq , 0) = phiScaHCurl(ishapeind , 0) * datavec[ hcurlmeshindex ].fNormalVec(0,ivecind);
        phiVecHDiv(iq , 1) = phiScaHCurl(ishapeind , 0) * datavec[ hcurlmeshindex ].fNormalVec(1,ivecind);
        phiVecHDiv(iq , 2) = phiScaHCurl(ishapeind , 0) * datavec[ hcurlmeshindex ].fNormalVec(2,ivecind);
    }
    
    /*********************CALCULATE NORMAL VECTOR****************************/
    TPZManVector<REAL,3> ax1(3),ax2(3), elNormal(3);
    for (int i=0; i<3; i++) {
        ax1[i] = datavec[ hcurlmeshindex ].axes(0,i);
        ax2[i] = datavec[ hcurlmeshindex ].axes(1,i);
    }
    Cross(ax1, ax2, elNormal);
    
    /*********************CREATE HCURL FUNCTIONS*****************************/
    TPZFNMatrix< 36 , REAL > phiHCurl(phrq , 3 , 0.);
    RotateForHCurl(elNormal , phiVecHDiv , phiHCurl);
    /*********************COMPUTE CURL****************************/
    TPZFMatrix<REAL> &dphiQdaxes = datavec[ hcurlmeshindex ].dphix;
    TPZFNMatrix<3,REAL> dphiQ;
    TPZAxesTools<REAL>::Axes2XYZ(dphiQdaxes, dphiQ, datavec[ hcurlmeshindex ].axes);
    TPZFNMatrix<3,REAL> gradPhiForHCurl(phrq , 3 , 0.);
    TPZFNMatrix<3,REAL> ivecHCurl(phrq , 3 , 0.);
    TPZManVector<REAL,3> iVecHDiv(3,0.), ivecForCurl(3,0.);
    for (int iPhi = 0; iPhi < phrq; iPhi++) {
        int ivecind = datavec[ hcurlmeshindex ].fVecShapeIndex[iPhi].first;
        int ishapeind = datavec[ hcurlmeshindex ].fVecShapeIndex[iPhi].second;
        iVecHDiv[0] = datavec[ hcurlmeshindex ].fNormalVec(0,ivecind);
        iVecHDiv[1] = datavec[ hcurlmeshindex ].fNormalVec(1,ivecind);
        iVecHDiv[2] = datavec[ hcurlmeshindex ].fNormalVec(2,ivecind);
        Cross(elNormal, iVecHDiv, ivecForCurl);
        for (int i = 0; i<dphiQ.Rows(); i++) {
            gradPhiForHCurl(iPhi,i) = dphiQ(i,ishapeind);
            ivecHCurl(iPhi,i) = ivecForCurl[i];
        }
    }
    TPZFNMatrix<40,REAL> curlPhi;
    ComputeCurl(gradPhiForHCurl, ivecHCurl, curlPhi);

    const REAL k0 = fScaleFactor * 2*M_PI/fLambda;
    /*****************ACTUAL COMPUTATION OF CONTRIBUTION****************/

    const int nHCurlFunctions  = phiHCurl.Rows();
    const int nH1Functions  = phiH1.Rows();
    const int firstH1 = h1meshindex * nHCurlFunctions;
    const int firstHCurl = hcurlmeshindex * nH1Functions;

    for (int iVec = 0; iVec < nHCurlFunctions; iVec++) {
        for (int jVec = 0; jVec < nHCurlFunctions; jVec++) {
            STATE curlIzdotCurlJz = 0.;
            curlIzdotCurlJz += curlPhi(iVec, 2) * curlPhi(jVec , 2);
            STATE phiIdotPhiJx = phiHCurl(iVec , 0) * phiHCurl(jVec , 0);
            STATE phiIdotPhiJy = phiHCurl(iVec , 1) * phiHCurl(jVec , 1);

            STATE stiffAtt = 0.;
            stiffAtt = (1./uzz) * curlIzdotCurlJz;
            stiffAtt -= k0 * k0 * exx * phiIdotPhiJx;
            stiffAtt -= k0 * k0 * eyy * phiIdotPhiJy;
            STATE stiffBtt = 0.;
            stiffBtt += (1./uyy) * phiIdotPhiJx;
            stiffBtt += (1./uxx) * phiIdotPhiJy;
            if (this->fAssembling == A) {
                ek( firstHCurl + iVec , firstHCurl + jVec ) += stiffAtt * weight ;
            }
            else if (this->fAssembling == B){
                ek( firstHCurl + iVec , firstHCurl + jVec ) += stiffBtt * weight ;
            }
            else{
                DebugStop();
            }

        }
        for (int jSca = 0; jSca < nH1Functions; jSca++) {
            STATE phiVecDotGradPhiScax = phiHCurl(iVec , 0) * gradPhiH1(jSca , 0);
            STATE phiVecDotGradPhiScay = phiHCurl(iVec , 1) * gradPhiH1(jSca , 1);

            STATE stiffBzt = 0.;
            stiffBzt += (1./uyy) * phiVecDotGradPhiScax;
            stiffBzt += (1./uxx) * phiVecDotGradPhiScay;
            if (this->fAssembling == A) {
                ek( firstHCurl + iVec , firstH1 + jSca ) += 0. ;
            }
            else if (this->fAssembling == B){
                ek( firstHCurl + iVec , firstH1 + jSca ) += stiffBzt * weight ;
            }
            else{
                DebugStop();
            }
        }
    }
    for (int iSca = 0; iSca < nH1Functions; iSca++) {
        for (int jVec = 0; jVec < nHCurlFunctions; jVec++) {
            STATE phiVecDotGradPhiScax = phiHCurl(jVec , 0) * gradPhiH1(iSca , 0);
            STATE phiVecDotGradPhiScay = phiHCurl(jVec , 1) * gradPhiH1(iSca , 1);

            STATE stiffBtz = 0.;
            stiffBtz += (1./uyy) * phiVecDotGradPhiScax;
            stiffBtz += (1./uxx) * phiVecDotGradPhiScay;
            if (this->fAssembling == A) {
                ek( firstH1 + iSca , firstHCurl +  jVec) += 0. ;
            }
            else if (this->fAssembling == B){
                ek( firstH1 + iSca , firstHCurl +  jVec ) += stiffBtz * weight ;
            }
            else{
                DebugStop();
            }
        }
        for (int jSca = 0; jSca < nH1Functions; jSca++) {
            STATE gradPhiScaDotGradPhiScax = gradPhiH1(iSca , 0) * gradPhiH1(jSca , 0);
            STATE gradPhiScaDotGradPhiScay = gradPhiH1(iSca , 1) * gradPhiH1(jSca , 1);

            STATE stiffBzz = 0.;
            stiffBzz +=  (1./uyy) * gradPhiScaDotGradPhiScax;
            stiffBzz +=  (1./uxx) * gradPhiScaDotGradPhiScay;
            stiffBzz -=  k0 * k0 * ezz * phiH1( iSca , 0 ) * phiH1( jSca , 0 );

            if (this->fAssembling == A) {
                ek( firstH1 + iSca , firstH1 + jSca) += 0. ;
            }
            else if (this->fAssembling == B){
                ek( firstH1 + iSca , firstH1 + jSca) += stiffBzz * weight ;
            }
            else{
                DebugStop();
            }
        }
    }
}



void TPZMatWaveguidePmlHDiv::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout)
{

    TPZVec<STATE> et(3,0.);
    TPZVec<STATE> ez(1,0.);
    TPZManVector<STATE,3> ax1(3),ax2(3), normal(3);
    for (int i=0; i<3; i++) {
        ax1[i] = datavec[hcurlmeshindex].axes(0,i);
        ax2[i] = datavec[hcurlmeshindex].axes(1,i);
    }
    //ROTATE FOR HCURL
    Cross(ax1, ax2, normal);
    
    Cross(normal,datavec[ hcurlmeshindex ].sol[0], et);
    ez = datavec[ h1meshindex ].sol[0];
    switch (var) {
        case 0:{//et
            Solout = et;
            break;
        }
        case 1:{//ez
            Solout = ez;
            break;
        }

        case 2:{//material
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
            Solout[0] = fEr * sx;
            Solout[1] = fEr * sy;
            break;
        }
        case 3:{//pOrder
            Solout.Resize(1);
            Solout[0] = datavec[h1meshindex].p;
            break;
        }
        case 4:{//pOrder
            Solout.Resize(1);
            Solout[0] = datavec[hcurlmeshindex].p;
            break;
        }
        default:
            DebugStop();
            break;
    }
}