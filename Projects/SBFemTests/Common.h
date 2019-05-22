#ifndef COMMONHPP
#define COMMONHPP

#include "pzanalysis.h"
#include "pzcmesh.h"
#include "TPZAnalyticSolution.h"

#ifdef _AUTODIFF
extern TElasticity2DAnalytic ElastExact;

extern TLaplaceExampleTimeDependent TimeLaplaceExact;

extern TLaplaceExample1 LaplaceExact;
#endif

//    Setup the system of equations and invert
void SolveSist(TPZAnalysis *an, TPZCompMesh *fCmesh);

/// insert material objects in the computational mesh
void InsertMaterialObjects(TPZCompMesh *cmesh, bool scalarproblem, bool applyexact);

void InsertMaterialObjects2(TPZCompMesh *cmesh, bool scalarproblem, bool applyexact);

/// Build a square mesh with boundary conditions
TPZCompMesh *SetupSquareMesh(int nelx, int nrefskeleton, int porder, bool elasticityproblem, bool applyexact);

/// Build a square mesh with boundary conditions
TPZCompMesh *SetupSquareH1Mesh(int nelx, int porder, bool elasticityproblem, bool applyexact);

/// Build a square mesh with boundary conditions
TPZCompMesh *SetupCrackedOneElement(int nrefskeleton, int porder, bool applyexact, bool elastic);

enum MMATID {Enomat, Emat1, Emat2, Emat3, Emat4, Ebc1, Ebc2, Ebc3, Ebc4, EBCPoint1, EBCPoint2, Ewrap, ESkeleton, EInterfaceMat1, EInterfaceMat2, EGroup};

//void Constant(const TPZVec<REAL> &x, TPZVec<STATE> &val);

/// Function defining the Harmonic solution at the left of the domain
void HarmonicNeumannLeft(const TPZVec<REAL> &x, TPZVec<STATE> &val);

/// Function defining the Harmonic solution at the right of the domain
void HarmonicNeumannRight(const TPZVec<REAL> &x, TPZVec<STATE> &val);

/// Function defining the exact harmonic solution
void Harmonic_exact(const TPZVec<REAL> &xv, TPZVec<STATE> &val, TPZFMatrix<STATE> &deriv);

#ifdef _AUTODIFF
/// Function defining the exact elasticity solution
inline void Elasticity_exact(const TPZVec<REAL> &xv, TPZVec<STATE> &val, TPZFMatrix<STATE> &deriv)
{
    ElastExact.Solution(xv, val, deriv);
}

inline void Laplace_exact(const TPZVec<REAL> &xv, TPZVec<STATE> &val, TPZFMatrix<STATE> &deriv)
{
    LaplaceExact.Solution(xv, val, deriv);
}

inline void TimeLaplace_exact(const TPZVec<REAL> &xv, TPZVec<STATE> &val, TPZFMatrix<STATE> &deriv)
{
    TimeLaplaceExact.Solution(xv, val, deriv);
}
#endif

/// Read a JSon File and generate a computational mesh
TPZCompMesh *ReadJSonFile(const std::string &filename, int numrefskeleton, int pOrder, REAL contrast);

/// Create a one element mesh going from angle = 0 to angle
TPZCompMesh *SetupOneArc(int numrefskeleton, int porder, REAL angle);

/// Verify if the values of the shapefunctions corresponds to the value of ComputeSolution for all SBFemVolumeElements
void VerifyShapeFunctionIntegrity(TPZCompMesh *cmesh);

void Constant_Laplacian_Ex1(const TPZVec<REAL> &x, TPZVec<STATE> &val);

void ExactSol_Laplacian_Ex1(const TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix<STATE> &deriv);

void Constant_Laplacian_Ex2(const TPZVec<REAL> &x, TPZVec<STATE> &val);

void ExactSol_Laplacian_Ex2(const TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix<STATE> &deriv);

void Constant_Laplacian_Ex3(const TPZVec<REAL> &x, TPZVec<STATE> &val);

void ExactSol_Laplacian_Ex3(const TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix<STATE> &deriv);

void Constant_Laplacian_Ex4(const TPZVec<REAL> &x, TPZVec<STATE> &val);

void ExactSol_Laplacian_Ex4(const TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix<STATE> &deriv);
#endif

void forcefunction(const TPZVec<REAL> &co, TPZVec<STATE> &result);

void solexact(const TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix<STATE> &deriv);


void WriteCSVFile(std::stringstream sout, TPZVec<REAL> &errors);

void SetExample(int example);

void SetInterpolation(TPZCompMesh * cmesh, REAL h);
