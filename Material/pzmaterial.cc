//METHODS DEFINITION FOR CLASS TPZMaterial

#include "pzmaterial.h"
#include "pzerror.h"
#include "pzvec.h"
#include "pzbndcond.h"
#include "pzreal.h"
#include "pzadmchunk.h"

REAL TPZMaterial::gBigNumber = 1.e12;

TPZMaterial::TPZMaterial(int id) {
   fId = id;
   fForcingFunction = 0;
}

TPZMaterial::TPZMaterial(TPZMaterial &material) {
   fId = material.fId;
   fForcingFunction = material.fForcingFunction;
}

void TPZMaterial::Print(ostream & out) {
   out << endl << "Material Id = " << fId << endl;
}

int TPZMaterial::VariableIndex(char *name) {
   if(!strcmp(name,"state")) return 0;
   if(!strcmp(name,"POrder")) return 99;
   if(!strcmp(name,"Error")) return 100;
   if(!strcmp(name,"TrueError")) return 101;
   if(!strcmp(name,"EffectivityIndex")) return 102;


   return -1;
}

int TPZMaterial::NSolutionVariables(int index) {
   if(index == 0) return NStateVariables();
   if(index == 99) return 1;
   if(index == 100) return 1;
   if(index == 101) return 1;
   if(index == 102) return 1;
   PZError << "TPZMaterial::NSolutionVariables called index = " << index << "\n";
   return 0;
}

void TPZMaterial::Solution(TPZVec<REAL> &Sol,TPZFMatrix &/*DSol*/,TPZFMatrix &/*axes*/,int var,
			   TPZVec<REAL> &Solout){
   if(var == 0) Solout = Sol;
   else if(var == 99 || var == 100 || var == 101 || var == 102) {
      //  	PZError << "TPZMaterial var = "<< var << " the element should treat this case\n";
      Solout[0] = Sol[0]; // = 0.;
   } else Solout.Resize(0);
}

TPZBndCond *TPZMaterial::CreateBC(int id, int typ, TPZFMatrix &val1, TPZFMatrix &val2) {
   return new TPZBndCond(this,id,typ,val1,val2);
}

void TPZMaterial::SetData(istream &data) {
   PZError << "TPZMaterial::SetData is called.\n";
   data >> fId;
}

TPZMaterial *TPZMaterial::NewMaterial() {
   PZError << "TPZMaterial::NewMaterial is called.\n";
   return 0;
}
void TPZMaterial::Contribute(TPZVec<REAL> &x,TPZFMatrix &jacinv, TPZVec<REAL> &sol,TPZFMatrix &dsol,REAL weight,TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &dphi,TPZFMatrix &ef){
   TPZFMatrix ek(ef.Rows(),ef.Rows(),0.);
   Contribute(x,jacinv,sol,dsol,weight,axes,phi,dphi,ek,ef);
}

void TPZMaterial::Clone(TPZAdmChunkVector<TPZMaterial *> &matvec) {
   int matid = Id();
   int nmat = matvec.NElements();
   int m;
   for(m=0; m<nmat; m++) {
      TPZMaterial *mat = matvec[m];
      if(!mat) continue;
      if(mat->Id() == matid) return;
   }
   int vecpos = matvec.AllocateNewElement();
   TPZMaterial *newmat = NewMaterial();
   matvec[vecpos] = newmat;
}
