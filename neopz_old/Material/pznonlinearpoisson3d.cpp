// -*- c++ -*- 

//$Id: pznonlinearpoisson3d.cpp,v 1.9 2007-12-05 14:15:30 tiago Exp $

#include "pznonlinearpoisson3d.h"
#include "pzbndcond.h"

using namespace std;

TPZNonLinearPoisson3d::TPZNonLinearPoisson3d(int nummat, int dim):TPZMatPoisson3dReferred(nummat, dim){
  this->fIsReferred = true;
  this->SetNoStabilizationTerm();
}

TPZNonLinearPoisson3d::TPZNonLinearPoisson3d(const TPZNonLinearPoisson3d &cp):TPZMatPoisson3dReferred(cp){
  this->fIsReferred = cp.fIsReferred;
  this->fStabilizationType = cp.fStabilizationType;
}

TPZNonLinearPoisson3d::~TPZNonLinearPoisson3d(){

}

void TPZNonLinearPoisson3d::SetSUPGStab(REAL sd){
  this->fStabilizationType = ESUPG;
  this->fSD = sd;
}
  
void TPZNonLinearPoisson3d::SetGradientStab(REAL sd){
  this->fStabilizationType = EGradient;
  this->fSD = sd;
}
  
void TPZNonLinearPoisson3d::SetNoStabilizationTerm(){
  this->fStabilizationType = ENoStabilization;
  this->fSD = 0.;
}

void TPZNonLinearPoisson3d::Contribute(TPZMaterialData &data,
                                       REAL weight,
                                       TPZFMatrix &ek,
                                       TPZFMatrix &ef){

  TPZFMatrix &dphi = data.dphix;
  TPZFMatrix &phi = data.phi;
  TPZManVector<REAL,3> &x = data.x;
  TPZVec<REAL> &sol=data.sol;
  TPZFMatrix &dsol=data.dsol;
  TPZFMatrix &axes=data.axes;
  TPZFMatrix &jacinv=data.jacinv;

  if (this->IsReferred()){
    this->SetConvectionTerm(dsol, axes);
  }

  int phr = phi.Rows();

  if(fForcingFunction) {      
    TPZManVector<REAL> res(1);
    fForcingFunction(x,res);  
    fXf = res[0];
  }
  REAL delx = 0.;
  REAL ConvDirAx[3] = {0.};
  if(fC != 0.0) {
    int di,dj;
    delx = 0.;
    for(di=0; di<fDim; di++) {
      for(dj=0; dj<fDim; dj++) {
        delx = (delx<fabs(jacinv(di,dj))) ? fabs(jacinv(di,dj)) : delx;
      }
    }
    delx = 2./delx;
      
    
    switch(fDim) {
      case 1:
        ConvDirAx[0] = axes(0,0)*fConvDir[0]+axes(0,1)*fConvDir[1]+axes(0,2)*fConvDir[2];
        break;
      case 2:
        ConvDirAx[0] = axes(0,0)*fConvDir[0]+axes(0,1)*fConvDir[1]+axes(0,2)*fConvDir[2];
        ConvDirAx[1] = axes(1,0)*fConvDir[0]+axes(1,1)*fConvDir[1]+axes(1,2)*fConvDir[2];
        break;
      case 3:
        ConvDirAx[0] = axes(0,0)*fConvDir[0]+axes(0,1)*fConvDir[1]+axes(0,2)*fConvDir[2];
        ConvDirAx[1] = axes(1,0)*fConvDir[0]+axes(1,1)*fConvDir[1]+axes(1,2)*fConvDir[2];
        ConvDirAx[2] = axes(2,0)*fConvDir[0]+axes(2,1)*fConvDir[1]+axes(2,2)*fConvDir[2];
        break;
      default:
        PZError << "TPZMatPoisson3d::Contribute dimension error " << fDim << endl;
    }
  }
 
  for( int in = 0; in < phr; in++ ) {
    int kd;
    REAL dphiic = 0;
    for(kd = 0; kd<fDim; kd++) dphiic += ConvDirAx[kd]*dphi(kd,in);
    ef(in, 0) += - weight * ( fXf*phi(in,0) + 0.5*fSD*delx*fC*dphiic*fXf );
    for(kd = 0; kd < fDim; kd++){
      ef(in, 0) += -1. * weight * ( +fK * ( dphi(kd,in) * dsol(kd,0) )
                                    -fC * ( ConvDirAx[kd]* dphi(kd,in) * sol[0] )  );
    }///kd

    for( int jn = 0; jn < phr; jn++ ) {
      for(kd=0; kd<fDim; kd++) {
        ek(in,jn) += weight * (
          +fK * ( dphi(kd,in) * dphi(kd,jn) ) 
          -fC * ( ConvDirAx[kd]* dphi(kd,in) * phi(jn) ) );
      }
    }
  }///in
    
  if (fStabilizationType == ESUPG){
  
    for( int in = 0; in < phr; in++ ) {
      int kd;
      REAL dphiic = 0;
      for(kd = 0; kd<fDim; kd++) dphiic += ConvDirAx[kd]*dphi(kd,in);
      ef(in, 0) += - weight * ( + 0.5*fSD*delx*fC*dphiic*fXf );
      for(kd = 0; kd < fDim; kd++){
        ef(in, 0) += -1. * weight * ( +0.5 * fSD * delx * fC * dphiic * dsol(kd,0) * ConvDirAx[kd] );
      }///kd
  
      for( int jn = 0; jn < phr; jn++ ) {
        for(kd=0; kd<fDim; kd++) {
          ek(in,jn) += weight * (            
            +0.5 * fSD * delx * fC * dphiic * dphi(kd,jn)* ConvDirAx[kd]
          );
        }
      }
    }///in
    
  }///SUPG
    
  if (fStabilizationType == EGradient){
    
    ///computing norm of solution gradient
    REAL dsolNorm = 0.;
    for(int d = 0; d < fDim; d++) dsolNorm += dsol(d,0)*dsol(d,0);
    dsolNorm = sqrt(dsolNorm);
    if (dsolNorm < 1e-16) dsolNorm = 1.;
    
    ///loop over i shape functions
    int kd;
    for( int in = 0; in < phr; in++ ){
    
      ///computing gradV.gradU/Norm(gradU)
      REAL dphiic = 0.;
      for(kd = 0; kd<fDim; kd++) dphiic += dsol(kd,0) * dphi(kd,in) / dsolNorm;

      ef(in, 0) += - weight * ( 0.5*fSD*delx* /*dsolNorm**/ dphiic* fXf );
      
      double aux = 0.;
      for(kd = 0; kd < fDim; kd++){
        aux += dphiic * dsol(kd,0)*(fC*ConvDirAx[kd]);
      }///kd
      ef(in,0) += -1.* ( +0.5 * fSD * delx * aux * weight );
  
      for( int jn = 0; jn < phr; jn++ ) {
        REAL DdphiicDalpha = 0.;
        for(kd = 0; kd<fDim; kd++) DdphiicDalpha += dphi(kd,jn)*dphi(kd,in)/dsolNorm;
        double aux = 0.;
        for(kd=0; kd<fDim; kd++) {
          aux += (fC*ConvDirAx[kd]) * ( dphiic*dphi(kd,jn) + dsol(kd,0)*DdphiicDalpha );
        }///kd
        ek(in,jn) += +0.5 * fSD * delx * aux * weight;        
      }///jn
      
    }///in
    
  }///EGradiente

  if (this->fC == 0.){
    if ( !ek.VerifySymmetry() ) cout << __PRETTY_FUNCTION__ << "\nMATRIZ NAO SIMETRICA" << endl;
  }

}

void TPZNonLinearPoisson3d::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix &ef){
  TPZFMatrix &dphi = data.dphix;
  TPZFMatrix &phi = data.phi;
  TPZManVector<REAL,3> &x = data.x;
  TPZVec<REAL> &sol=data.sol;
  TPZFMatrix &dsol=data.dsol;
  TPZFMatrix &axes=data.axes;
  TPZFMatrix &jacinv=data.jacinv;

  if (this->IsReferred()){
    this->SetConvectionTerm(dsol, axes);
  }

  int phr = phi.Rows();

  if(fForcingFunction) {      
    TPZManVector<REAL> res(1);
    fForcingFunction(x,res);  
    fXf = res[0];
  }
  REAL delx = 0.;
  REAL ConvDirAx[3] = {0.};
  if(fC != 0.0) {
    int di,dj;
    delx = 0.;
    for(di=0; di<fDim; di++) {
      for(dj=0; dj<fDim; dj++) {
        delx = (delx<fabs(jacinv(di,dj))) ? fabs(jacinv(di,dj)) : delx;
      }
    }
    delx = 2./delx;
      
    
    switch(fDim) {
      case 1:
        ConvDirAx[0] = axes(0,0)*fConvDir[0]+axes(0,1)*fConvDir[1]+axes(0,2)*fConvDir[2];
        break;
      case 2:
        ConvDirAx[0] = axes(0,0)*fConvDir[0]+axes(0,1)*fConvDir[1]+axes(0,2)*fConvDir[2];
        ConvDirAx[1] = axes(1,0)*fConvDir[0]+axes(1,1)*fConvDir[1]+axes(1,2)*fConvDir[2];
        break;
      case 3:
        ConvDirAx[0] = axes(0,0)*fConvDir[0]+axes(0,1)*fConvDir[1]+axes(0,2)*fConvDir[2];
        ConvDirAx[1] = axes(1,0)*fConvDir[0]+axes(1,1)*fConvDir[1]+axes(1,2)*fConvDir[2];
        ConvDirAx[2] = axes(2,0)*fConvDir[0]+axes(2,1)*fConvDir[1]+axes(2,2)*fConvDir[2];
        break;
      default:
        PZError << "TPZMatPoisson3d::Contribute dimension error " << fDim << endl;
    }
  }
 
  for( int in = 0; in < phr; in++ ) {
    int kd;
    REAL dphiic = 0;
    for(kd = 0; kd<fDim; kd++) dphiic += ConvDirAx[kd]*dphi(kd,in);
    ef(in, 0) += - weight * ( fXf*phi(in,0) + 0.5*fSD*delx*fC*dphiic*fXf );
    for(kd = 0; kd < fDim; kd++){
      ef(in, 0) += -1. * weight * ( +fK * ( dphi(kd,in) * dsol(kd,0) )
                                    -fC * ( ConvDirAx[kd]* dphi(kd,in) * sol[0] )  );
    }///kd

  }///in
    
  if (fStabilizationType == ESUPG){
  
    for( int in = 0; in < phr; in++ ) {
      int kd;
      REAL dphiic = 0;
      for(kd = 0; kd<fDim; kd++) dphiic += ConvDirAx[kd]*dphi(kd,in);
      ef(in, 0) += - weight * ( + 0.5*fSD*delx*fC*dphiic*fXf );
      for(kd = 0; kd < fDim; kd++){
        ef(in, 0) += -1. * weight * ( +0.5 * fSD * delx * fC * dphiic * dsol(kd,0) * ConvDirAx[kd] );
      }///kd
  
    }///in
    
  }///SUPG
    
  if (fStabilizationType == EGradient){
    
    ///computing norm of solution gradient
    REAL dsolNorm = 0.;
    for(int d = 0; d < fDim; d++) dsolNorm += dsol(d,0)*dsol(d,0);
    dsolNorm = sqrt(dsolNorm);
    if (dsolNorm < 1e-16) dsolNorm = 1.;
    
    ///loop over i shape functions
    int kd;
    for( int in = 0; in < phr; in++ ){
    
      ///computing gradV.gradU/Norm(gradU)
      REAL dphiic = 0.;
      for(kd = 0; kd<fDim; kd++) dphiic += dsol(kd,0) * dphi(kd,in) / dsolNorm;

      ef(in, 0) += - weight * ( 0.5*fSD*delx* /*dsolNorm**/ dphiic* fXf );
      
      double aux = 0.;
      for(kd = 0; kd < fDim; kd++){
        aux += dphiic * dsol(kd,0)*(fC*ConvDirAx[kd]);
      }///kd
      ef(in,0) += -1.* ( +0.5 * fSD * delx * aux * weight );

    }///in
    
  }///EGradiente

}///void
               
void TPZNonLinearPoisson3d::ContributeBC(TPZMaterialData &data,
                                         REAL weight,
                                         TPZFMatrix &ek,
                                         TPZFMatrix &ef,
                                         TPZBndCond &bc){
  TPZFMatrix &phi = data.phi;
  TPZVec<REAL> &sol=data.sol;
  TPZFMatrix &axes=data.axes;

  int phr = phi.Rows();
  short in,jn;
  REAL v2[1];
  v2[0] = bc.Val2()(0,0);

  switch (bc.Type()) {
  case 0 : {      // Dirichlet condition
    for(in = 0 ; in < phr; in++) {
      ef(in,0) += weight * ( gBigNumber * v2[0] * phi(in,0) - gBigNumber * phi(in,0) * sol[0] );
      for (jn = 0 ; jn < phr; jn++) {
        ek(in,jn) += gBigNumber * phi(in,0) * phi(jn,0) * weight;
      }
    }
  }
  break;
    
  case 1 : {      // Neumann condition
    for(in = 0 ; in < phi.Rows(); in++) {
      ef(in,0) += v2[0] * phi(in,0) * weight;
    }
  }
  break;
    
  case 2 :{    // condicao mista
    cout << __PRETTY_FUNCTION__ << " at line " << __LINE__ << " not implemented\n";
  }
  break;
    
  case 3 : { // outflow condition
  
    if (this->IsReferred()){
      PZError << "Error at " << __PRETTY_FUNCTION__
              << " - the outflow boundary condition can not be implemented for referred elements derived from TPZInterpolatedElement\n";
    }
  
    int id, il, jl;
    REAL normal[3];
    if (fDim == 1) PZError << __PRETTY_FUNCTION__ << " - ERROR! The normal vector is not available for 1D TPZInterpolatedElement\n";
    if (fDim == 2){
      normal[0] = axes(0,1);
      normal[1] = axes(1,1);
    }
    if (fDim == 3){
      normal[0] = axes(0,2);
      normal[1] = axes(1,2);
      normal[2] = axes(2,2);
    }
    REAL ConvNormal = 0.;    
    for(id=0; id<fDim; id++) ConvNormal += fC*fConvDir[id]*normal[id];  
    if(ConvNormal > 0.) {
      for(il=0; il<phr; il++) {
        for(jl=0; jl<phr; jl++) {
          ek(il,jl) += weight * ConvNormal * phi(il)*phi(jl);
        }
        ef(il,0) += -1. * weight * ConvNormal * phi(il) * sol[0];
      }
    }
    else{
      if (ConvNormal < 0.) std::cout << "Boundary condition error: inflow detected in outflow boundary condition: ConvNormal = " << ConvNormal << "\n";
    }  
  }
  break;
    
  default :{
    PZError << __PRETTY_FUNCTION__ << " at line " << __LINE__ << " - Error! Wrong boundary condition type\n";
  }
  break;
  }
  
  if (this->IsSymetric()) {//only 1.e-3 because of bignumbers.
    if ( !ek.VerifySymmetry( 1.e-3 ) ) cout << __PRETTY_FUNCTION__ << "\nMATRIZ NAO SIMETRICA" << endl;
  }
}

void TPZNonLinearPoisson3d::ContributeInterface(TPZMaterialData &data,
                                                REAL weight,
                                                TPZFMatrix &ek,
                                                TPZFMatrix &ef){
TPZFMatrix &dphiL = data.dphixl;
TPZFMatrix &dphiR = data.dphixr;
TPZFMatrix &phiL = data.phil;
TPZFMatrix &phiR = data.phir;
TPZManVector<REAL,3> &normal = data.normal;
TPZVec<REAL> &solL=data.soll;
TPZVec<REAL> &solR=data.solr;
TPZFMatrix &dsolL=data.dsoll;
TPZFMatrix &dsolR=data.dsolr;

  if (this->IsReferred()){
    this->SetConvectionTermInterface(dsolL, dsolR);
  }
  
  const int nrowl = phiL.Rows();
  const int nrowr = phiR.Rows();
  
  ///Convection term
  REAL ConvNormal = 0.;
  for(int id=0; id<fDim; id++) ConvNormal += fC * fConvDir[id]*normal[id];
  if(ConvNormal > 0.) {
    for(int il=0; il<nrowl; il++) {
      ef(il, 0) += -1. * weight * ConvNormal * phiL(il) * solL[0];
      for(int jl=0; jl<nrowl; jl++) {
        ek(il,jl) += weight * ConvNormal * phiL(il)*phiL(jl);
      }
    }
    for(int ir=0; ir<nrowr; ir++) {
      ef(ir+nrowl,0) += -1. * (-1. * weight * ConvNormal * phiR(ir) * solL[0] );
      for(int jl=0; jl<nrowl; jl++) {
        ek(ir+nrowl,jl) -= weight * ConvNormal * phiR(ir) * phiL(jl);
      }
    }
  } else {
    for(int ir=0; ir<nrowr; ir++) {
      ef(ir+nrowl,0) += -1. * (-1. * weight * ConvNormal * phiR(ir) * solR[0] );
      for(int jr=0; jr<nrowr; jr++) {
        ek(ir+nrowl,jr+nrowl) -= weight * ConvNormal * phiR(ir) * phiR(jr);
      }
    }
    for(int il=0; il<nrowl; il++) {
      ef(il,0) += -1. * weight * ConvNormal * phiL(il) * solR[0];
      for(int jr=0; jr<nrowr; jr++) {
        ek(il,jr+nrowl) += weight * ConvNormal * phiL(il) * phiR(jr);
      }
    }
  }
  
  
  ///diffusion term
  REAL leftK, rightK;
  leftK  = this->fK;
  rightK = this->GetRightK();
  
  //Compute GradSol . normal
  REAL DSolLNormal = 0.;
  REAL DSolRNormal = 0.;
  for(int id=0; id<fDim; id++) {
    DSolLNormal += dsolL(id,0)*normal[id];
    DSolRNormal += dsolR(id,0)*normal[id];
  }//for

  // 1) phi_I_left, phi_J_left
  for(int il=0; il<nrowl; il++) {
    REAL dphiLinormal = 0.;
    for(int id=0; id<fDim; id++) {
      dphiLinormal += dphiL(id,il)*normal[id];
    }
    
    //ef = F - K u
    ef(il,0) += -1. * (weight * leftK * (this->fSymmetry * 0.5 * dphiLinormal*solL[0]-0.5*DSolLNormal*phiL(il,0)));
    
    for(int jl=0; jl<nrowl; jl++) {
      REAL dphiLjnormal = 0.;
      for(int id=0; id<fDim; id++) {
        dphiLjnormal += dphiL(id,jl)*normal[id];
      }
      ek(il,jl) += weight * leftK * (
        this->fSymmetry * 0.5*dphiLinormal*phiL(jl,0)-0.5*dphiLjnormal*phiL(il,0)
      );
    }
  }
  
  // 2) phi_I_right, phi_J_right
  for(int ir=0; ir<nrowr; ir++) {
    REAL dphiRinormal = 0.;
    for(int id=0; id<fDim; id++) {
      dphiRinormal += dphiR(id,ir)*normal[id];
    }

    //ef = F - K u
    ef(ir+nrowl,0) += -1. * weight * rightK * ( this->fSymmetry * (-0.5 * dphiRinormal * solR[0] ) + 0.5 * DSolRNormal * phiR(ir) );
    
    for(int jr=0; jr<nrowr; jr++) {
      REAL dphiRjnormal = 0.;
      for(int id=0; id<fDim; id++) {
        dphiRjnormal += dphiR(id,jr)*normal[id];
      }
      ek(ir+nrowl,jr+nrowl) += weight * rightK * (
        this->fSymmetry * (-0.5 * dphiRinormal * phiR(jr) ) + 0.5 * dphiRjnormal * phiR(ir)
      );
    }
  }
  
  // 3) phi_I_left, phi_J_right
  for(int il=0; il<nrowl; il++) {
    REAL dphiLinormal = 0.;
    for(int id=0; id<fDim; id++) {
      dphiLinormal += dphiL(id,il)*normal[id];
    }
    
    //ef = F - K u
    ef(il,0) += -1. * weight * ( this->fSymmetry * (-0.5 * dphiLinormal * leftK * solR[0] ) - 0.5 * DSolRNormal * rightK * phiL(il) );
    
    for(int jr=0; jr<nrowr; jr++) {
      REAL dphiRjnormal = 0.;
      for(int id=0; id<fDim; id++) {
        dphiRjnormal += dphiR(id,jr)*normal[id];
      }
      ek(il,jr+nrowl) += weight * (
        this->fSymmetry * (-0.5 * dphiLinormal * leftK * phiR(jr) ) - 0.5 * dphiRjnormal * rightK * phiL(il)
      );
    }
  }
  
  // 4) phi_I_right, phi_J_left
  for(int ir=0; ir<nrowr; ir++) {
    REAL dphiRinormal = 0.;
    for(int id=0; id<fDim; id++) {
      dphiRinormal += dphiR(id,ir)*normal[id];
    }
    
    //ef = F - K u
    ef(ir+nrowl,0) += -1. * weight * (this->fSymmetry * 0.5 * dphiRinormal * rightK * solL[0] + 0.5 * DSolLNormal * leftK * phiR(ir));
    
    for(int jl=0; jl<nrowl; jl++) {
      REAL dphiLjnormal = 0.;
      for(int id=0; id<fDim; id++) {
        dphiLjnormal += dphiL(id,jl)*normal[id];
      }
      ek(ir+nrowl,jl) += weight * (
        this->fSymmetry * 0.5 * dphiRinormal * rightK * phiL(jl) + 0.5 * dphiLjnormal * leftK * phiR(ir)
      );
    }
  }
  
  if (this->IsSymetric()){
    if ( !ek.VerifySymmetry() ) cout << __PRETTY_FUNCTION__ << "\nMATRIZ NAO SIMETRICA" << endl;
  }
}

void TPZNonLinearPoisson3d::ContributeBCInterface(TPZMaterialData &data,
                                                  REAL weight, 
                                                  TPZFMatrix &ek,
                                                  TPZFMatrix &ef,
                                                  TPZBndCond &bc) {
TPZFMatrix &dphiL = data.dphixl;
TPZFMatrix &phiL = data.phil;
TPZManVector<REAL,3> &normal = data.normal;
TPZVec<REAL> &solL=data.soll;
TPZFMatrix &dsolL=data.dsoll;

if (this->IsReferred()){
    this->SetConvectionTermInterface(dsolL, dsolL);
  }
  
  int il,jl,nrowl,id;
  nrowl = phiL.Rows();
  REAL ConvNormal = 0.;
  for(id=0; id<fDim; id++) ConvNormal += fC * fConvDir[id]*normal[id];
  
  //Compute GradSol . normal
  REAL DSolLNormal = 0.;
  for(id=0; id<fDim; id++) {
    DSolLNormal += dsolL(id,0)*normal[id];
  }//for
  
  switch(bc.Type()) {
  case 0: // DIRICHLET
    
    //Diffusion
    for(il=0; il<nrowl; il++) {
      REAL dphiLinormal = 0.;
      for(id=0; id<fDim; id++) {
        dphiLinormal += dphiL(id,il)*normal[id];
      }
      ef(il,0) += weight*fK*dphiLinormal*bc.Val2()(0,0) * this->fSymmetry;
      
      //ef = F - K u
      ef(il,0) += -1. * weight*fK*(this->fSymmetry * dphiLinormal * solL[0] - DSolLNormal * phiL(il,0));
      
      for(jl=0; jl<nrowl; jl++) {
        REAL dphiLjnormal = 0.;
        for(id=0; id<fDim; id++) {
          dphiLjnormal += dphiL(id,jl)*normal[id];
        }
        ek(il,jl) += weight*fK*(this->fSymmetry * dphiLinormal * phiL(jl,0) - dphiLjnormal * phiL(il,0));
      }
    }
    
    //Convection
    if(ConvNormal > 0.) {
      for(il=0; il<nrowl; il++) {
        ef(il,0) += -1. * weight * ConvNormal * phiL(il) * solL[0];
        for(jl=0; jl<nrowl; jl++) {
          ek(il,jl) += weight * ConvNormal * phiL(il)*phiL(jl);
        }
      }
    } else {
      for(il=0; il<nrowl; il++) {
        ef(il,0) -= weight * ConvNormal * bc.Val2()(0,0) * phiL(il);
      }
    }
      
    break;
  case 1: // Neumann
    for(il=0; il<nrowl; il++) {
      ef(il,0) += weight*phiL(il,0)*bc.Val2()(0,0);
    }
    break;
    
  case 3: // outflow condition
    if(ConvNormal > 0.) {
      for(il=0; il<nrowl; il++) {
        for(jl=0; jl<nrowl; jl++) {
          ek(il,jl) += weight * ConvNormal * phiL(il)*phiL(jl);
        }
        ef(il,0) += -1. * weight * ConvNormal * phiL(il) * solL[0];
      }
    }
    else {
      if (ConvNormal < 0.){
        std::cout << "Boundary condition error: inflow detected in outflow boundary condition: ConvNormal = " << ConvNormal << "\n";
      }
    }
    break;    
    
  default:
    PZError << __PRETTY_FUNCTION__ << " - Wrong boundary condition type\n";
    break;
  }
    if (this->IsSymetric()){
      if ( !ek.VerifySymmetry() ) cout << __PRETTY_FUNCTION__ << "\nMATRIZ NAO SIMETRICA" << endl;
    }
}
