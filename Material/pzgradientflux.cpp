//$Id: pzgradientflux.cpp,v 1.1 2009-08-28 19:43:44 fortiago Exp $

#include "pzgradientflux.h"

TPZGradientFlux::TPZGradientFlux(){

}

TPZGradientFlux::TPZGradientFlux(const TPZGradientFlux &cp){

}

TPZGradientFlux::~TPZGradientFlux(){

}

void TPZGradientFlux::ComputeFlux(TPZVec<REAL> &solL, TPZVec<REAL> &solR, const TPZVec<REAL> &normal, TPZVec<REAL> & F){
  const int nstate = 5;
  const int dim = 3;
  F.Resize(nstate*dim);
  F.Fill(0.);
  for(int i = 0; i < nstate; i++){
    for(int id = 0; id < dim; id++){
      F[i*dim+id] = 0.5*(solL[i]+solR[i])*normal[id];
    }///for id
  }///for i
}///void

void TPZGradientFlux::ApplyLimiter(TPZMaterialData &data){
  const int nstate = 5;
  const int dim = 3;
  TPZManVector<REAL,dim> dL(dim), dR(dim);
  for(int id = 0; id < dim; id++){
    dL[id] = data.x[id] - data.XLeftElCenter[id];
    dR[id] = data.x[id] - data.XRightElCenter[id];
  }///for id

  TPZManVector<REAL,dim> gradL(3), gradR(3);
  for(int is = 0; is < nstate; is++){
    for(int id = 0; id < dim; id++){
      gradL[id] = data.soll[ nstate + is*dim +id ];
      gradR[id] = data.solr[ nstate + is*dim +id ];
    }///for id
    this->ApplyVanAlbadaLimiter(data.soll[is], data.solr[is], gradL, gradR, data.normal, dL, dR);
//     this->ApplyMinModLimiter(data.soll[is], data.solr[is], gradL, gradR, data.normal, dL, dR);
  }///for is
}///void

void TPZGradientFlux::ApplyMinModLimiter(REAL &soll, REAL &solr,
                                         const TPZVec<REAL>& gradL, const TPZVec<REAL> &gradR,
                                         const TPZVec<REAL> &normal, 
                                         const TPZVec<REAL> &dL, const TPZVec<REAL> & dR){
  const int dim = 3;
  const double sL = this->Dot(gradL,normal);
  const double sR = this->Dot(gradR,normal);
  if((sL*sR) < 0.){
    return;
  }
  const double k = 0.5;
  if(fabs(sL) < fabs(sR)){
    soll = soll + k*this->Dot(gradL,dL);
    solr = solr + k*this->Dot(gradL,dR);
    return;
  }
  if(fabs(sL) > fabs(sR)){
    soll = soll + k*this->Dot(gradR,dL);
    solr = solr + k*this->Dot(gradR,dR);
    return;
  }
  ///if (sL == sR)
  soll = soll + k*this->Dot(gradL,dL);
  solr = solr + k*this->Dot(gradR,dR);
}///void

void TPZGradientFlux::ApplyVanAlbadaLimiter(REAL &soll, REAL &solr,
                                            const TPZVec<REAL>& gradL, const TPZVec<REAL> &gradR,
                                            const TPZVec<REAL> &normal, 
                                            const TPZVec<REAL> &dL, const TPZVec<REAL> & dR){
  const int dim = 3;
  const double sL = this->Dot(gradL,normal);
  const double sR = this->Dot(gradR,normal);
  if(fabs(sL) < 1e-12){
    return;
  }
  const double r = sR/sL;
  double phi = (r <= 0.) ? 0. : ((r*r+r)/(1.+r*r));
  phi *= 0.5;
  soll = soll + phi*this->Dot(gradL,dL);
  solr = solr + phi*this->Dot(gradR,dR);

}///void
