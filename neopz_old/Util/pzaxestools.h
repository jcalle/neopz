//$Id: pzaxestools.h,v 1.1 2008-02-12 10:35:31 tiago Exp $

// -*- c++ -*-

#ifndef AXESTOOLS
#define AXESTOOLS

#include "pzfmatrix.h"
#include "pzerror.h"

class TPZAxesTools{ 
  public:

  TPZAxesTools(){}

  ~TPZAxesTools(){}

  /** Verify whether parameter axes is an orthogonal normalized matrix. */
  static void VerifyAxes(const TPZFMatrix &axes){
#ifdef DEBUG
    const double tol = 1.e-7;
    bool check = true;
    for(int i = 0; i < axes.Rows(); i++){
      for(int j = 0; j < axes.Rows(); j++){
        if(i == j) continue;
        double sum = 0.;
        for(int k = 0; k < axes.Cols(); k++){
          sum += axes.GetVal(i,k)*axes.GetVal(j,k);
        }///k
        if(sum > tol){
          check = false;
        }///if
      }///for j
    }///for i    
    
    if(check == false){
      PZError << "\nError at " << __PRETTY_FUNCTION__;
      axes.Print(", axes = ", PZError);
      PZError <<  "\n";
    }
#endif
  }

  /** Makes the basis transformation from axes basis to euclidian basis.
    * @param axes must be an orthogonal normalized matrix. Axes vectors are written in rows.
    */
  static void Axes2XYZ(const TPZFMatrix &dudaxes, TPZFMatrix &dudx, const TPZFMatrix &axes){
    TPZAxesTools::VerifyAxes(axes);
    TPZFNMatrix<9> axesT;
    axes.Transpose(&axesT);
    dudx.Resize(axesT.Rows(), dudaxes.Cols());
    dudx.Zero();
    axesT.Multiply(dudaxes,dudx);
  }

  /** Makes the basis transformation from euclidian basis to axes basis 
    * @param axes must be an orthogonal normalized matrix. Axes vectors are written in rows.
    */
  static void XYZ2Axes(TPZFMatrix &dudaxes, const TPZFMatrix &dudx, const TPZFMatrix &axes){
    TPZAxesTools::VerifyAxes(axes);
    dudaxes.Resize(axes.Rows(), dudx.Cols());
    dudaxes.Zero();
    axes.Multiply(dudx,dudaxes);  
  }
  
  /** Test code
   */
  static int main(){
    TPZFMatrix axes(2,3,0.);
    axes(0,0) = sqrt(2.)/2.;
    axes(0,1) = 0.;
    axes(0,2) = sqrt(2.)/2.;
    axes(1,0) = 0.;
    axes(1,1) = 1.;
    axes(1,2) = 0.;
    TPZFMatrix dudaxes(2,1);
    dudaxes(0,0) = 0.3;
    dudaxes(1,0) = 0.7;
    TPZFMatrix dudx;
    dudaxes.Print("dudaxes=",std::cout);
    TPZAxesTools::Axes2XYZ(dudaxes,dudx,axes);
    dudx.Print("dudx=",std::cout);
    TPZAxesTools::XYZ2Axes(dudaxes,dudx,axes);
    dudaxes.Print("dudaxes=",std::cout);
    return 1;
}

};

#endif
