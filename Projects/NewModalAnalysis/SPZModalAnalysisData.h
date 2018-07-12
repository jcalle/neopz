//
// Created by Francisco Teixeira Orlandini on 2/8/18.
//

#ifndef PZ_SPZMODALANALYSISDATA_H
#define PZ_SPZMODALANALYSISDATA_H

#include <slepceps.h>
#include <pzreal.h>
#include <pzvec.h>
#include "parameter_handler.h"

struct SPZModalAnalysisData{
  enum PzCases{
    StepFiber = 1, RectangularWG = 2, HoleyFiber
  };
  enum boundtype{
    PEC = 0, PMC = 1
  };
  enum pmltype{
    xp=0,yp,xm,ym,xpyp,xmyp,xmym,xpym
  };
  struct SPZRectangularWGOpts{
    REAL wDomain;
    REAL hDomain;
    bool usingSymmetry;
    boundtype symmetryType;
  };
  struct SPZStepFiberOpts{
    REAL realRCore;
    REAL dPML;
    REAL boundDist;
    REAL hasPML;
    int nLayersPml;
  };

  struct SPZHoleyFiberOpts{
      REAL dPML;
      REAL boundDist;
      SPZModalAnalysisData::boundtype symmetryX;
      SPZModalAnalysisData::boundtype symmetryY;
      bool refineH;
      bool refineP;
  };
  struct SPZPhysicalOpts{
    bool isCutOff;
    REAL lambda;
    int nMaterials;
    TPZVec<STATE> urVec;
    TPZVec<STATE> erVec;
    TPZVec<REAL> freqVec;
    bool isLambda;
    REAL alphaMax;
    SPZRectangularWGOpts rectangularWgOpts;
    SPZStepFiberOpts stepFiberOpts;
    SPZHoleyFiberOpts holeyFiberOpts;
    TPZVec<pmltype> pmlTypeVec;
    TPZVec<boundtype> boundTypeVec;
  };
  struct SPZPzOpts{
// polynomial order of basis functions
      int pOrder;
// generate vtk for fields visualisation
      bool genVTK;
//whether to calculate error analysis
      bool l2error;
// export l2 error
      bool exportl2error;
// export eigen values
      bool exportEigen;
//number of NeoPZ threads
      int nThreads;
//prefix to be added to all exported files
      std::string prefix;
      //whether to calculate abs or re of the eigenvectors
      bool absVal;
      //vtk resolution
      long vtkRes;
      bool exportGMesh;
      bool exportCMesh;
      REAL scaleFactor;
      bool isTargetScaled;
      std::string meshFile;
      bool externGenMesh;
      int meshOrder;
      int pSteps;
      int hSteps;
      TPZVec<int> factorVec;
      bool scaleByk0;
      bool usingNeoPzMesh;
      PzCases pzCase;
  };
  struct SPZSolverOpts{
      EPSProblemType eps_prob_type;
      EPSType eps_type;
      bool eps_krylov_locking;
      PetscReal eps_krylov_restart;
      EPSConv eps_conv_test;
      bool eps_true_res;
      EPSWhich eps_which_eig;
      PetscScalar target;
      PetscReal eps_tol;
      PetscInt eps_max_its;
      PetscInt eps_nev;
      PetscInt eps_ncv;
      PetscInt eps_mpd;
      PetscInt eps_verbose;

      PCType st_precond;
      KSPType st_solver;
      PetscReal ksp_rtol;
      PetscReal ksp_atol;
      PetscReal ksp_dtol;
      PetscReal ksp_max_its;
      STType st_type;
  };

  SPZPhysicalOpts physicalOpts;
  SPZPzOpts pzOpts;
  SPZSolverOpts solverOpts;
};

#endif //PZ_SPZMODALANALYSISDATA_H
