/*T
Examples based on PETSc tutorials
*/

#include "tutorials.h"
#include <petscksp.h>
#include <pzerror.h>
#include "TPZPetScMatrix.h"

#include "TPZGenGrid2D.h"
#include "pzcmesh.h"
#include "TPZMatLaplacian.h"
#include "TPZAnalyticSolution.h"
#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "pzskylstrmatrix.h"
#include "pzbndcond.h"

#ifdef _AUTODIFF
TLaplaceExample1 LaplaceExact;
inline void Laplace_exact(const TPZVec<REAL> &xv, TPZVec<STATE> &val, TPZFMatrix<STATE> &deriv)
{
    LaplaceExact.Solution(xv, val, deriv);
}
#endif

enum MMATID {Enomat, Emat1, Ebc1, Ebc2, Ebc3, Ebc4};

TPZGeoMesh * SetupGeom(int &nelx);

void InsertMaterialObjects(TPZCompMesh *cmesh);

void SolveSist(TPZAnalysis * an, TPZCompMesh *Cmesh, int &numthreads);

int main(int argc,char **args)
{
#ifdef _AUTODIFF
    LaplaceExact.fExact = TLaplaceExample1::ECosCos;
#endif

    // // Just testing the class:
    // static char help[] = "Testing TPZPetscMatrix";
    // int ierr = PetscInitialize(&argc,&args,(char*)0,help); if (ierr) return 0;

    int nelx = 2, porder = 3, numthreads = 4;
    TPZGeoMesh * gmesh = SetupGeom(nelx);

    TPZCompMesh *FEM = new TPZCompMesh(gmesh);
    FEM->SetDefaultOrder(porder);
    InsertMaterialObjects(FEM);
    FEM->AutoBuild();
    int64_t neq = FEM->NEquations();
    std::cout << "neq = " << neq << "\n";

    bool mustOptimizeBandwidth = true;
    TPZAnalysis * an = new TPZAnalysis(FEM,mustOptimizeBandwidth);
    SolveSist(an,FEM,numthreads);

#ifdef _AUTODIFF
    an->SetExact(Laplace_exact);
#endif
    TPZStack<std::string> vecnames,scalnames;
    scalnames.Push("State");
    an->DefineGraphMesh(2, scalnames, vecnames, "../RegularSolutionfilename.vtk");
    an->PostProcess(3);
    
    std::cout << "Compute errors\n";
    TPZManVector<REAL,10> errors(3,0.);
    an->SetThreadsForError(numthreads);
    an->PostProcessError(errors,false);

    delete FEM;
    delete an;

    // PetscFinalize();

    // Testing the systems:
    // int ierr = tridiaglinsystem(argc, args);CHKERRQ(ierr);
    // ierr = tridiaglinsystemMUMPS(argc, args);CHKERRQ(ierr);
    // ierr = parallellinsystem(argc, args);CHKERRQ(ierr);

}

void InsertMaterialObjects(TPZCompMesh *cmesh)
{
    TPZMaterial *material;
    TPZMatLaplacian *matloc = new TPZMatLaplacian(Emat1);
    matloc->SetDimension(2);
    matloc->SetSymmetric();
    material = matloc;

    TPZFMatrix<STATE> val1(1, 1, 0.), val2(1, 1, 0.);
    TPZMaterial *BCond1 = material->CreateBC(material, Ebc1, 0, val1, val2);
    TPZMaterial *BCond2 = material->CreateBC(material, Ebc2, 0, val1, val2);
    TPZMaterial *BCond3 = material->CreateBC(material, Ebc3, 0, val1, val2);
    TPZMaterial *BCond4 = material->CreateBC(material, Ebc4, 0, val1, val2);
#ifdef _AUTODIFF
    BCond1->SetForcingFunction(LaplaceExact.TensorFunction());
    BCond2->SetForcingFunction(LaplaceExact.TensorFunction());
    BCond3->SetForcingFunction(LaplaceExact.TensorFunction());
    BCond4->SetForcingFunction(LaplaceExact.TensorFunction());
#endif

    cmesh->InsertMaterialObject(material);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
}

TPZGeoMesh * SetupGeom(int &nelx)
{
    TPZManVector<REAL, 4> x0(3, -1.), x1(3, 1.);
    x0[0] = -1, x0[1] = -1, x0[2] = 0.;
    x1[0] = 1, x1[1] = 1, x1[2] = 0.;

    TPZManVector<int, 4> nx(2, nelx);
    TPZGenGrid2D gengrid(nx, x0, x1);
    gengrid.SetElementType(MMeshType::EQuadrilateral);
    TPZGeoMesh * gmesh = new TPZGeoMesh();
    gengrid.Read(gmesh, Emat1);
    gengrid.SetBC(gmesh, 4, Ebc1);
    gengrid.SetBC(gmesh, 5, Ebc2);
    gengrid.SetBC(gmesh, 6, Ebc3);
    gengrid.SetBC(gmesh, 7, Ebc4);
    gmesh->BuildConnectivity();

    return gmesh;
}

void SolveSist(TPZAnalysis * an, TPZCompMesh *Cmesh, int &numthreads)
{
    TPZPetScMatrix<REAL> strmat(Cmesh);
    // TPZSkylineStructMatrix strmat(Cmesh);
    // strmat.SetNumThreads(numthreads);
    an->SetStructuralMatrix(strmat);

    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky);
    an->SetSolver(step);

    an->Run();
}

void ClassTests()
{
    // TPZPetScMatrix<REAL> mat();
    // TPZPetScMatrix<REAL> mat2(2,2);
    TPZPetScMatrix<REAL> mat1(2,2,4);
    mat1.Print("m1=", std::cout, EMathematicaInput);
    TPZPetScMatrix<REAL> mat2(2,2,5);
    mat2.Print("m2=", std::cout, EMathematicaInput);
    mat1.Zero();
    mat1(1,1) = 3;
    mat1.Print("m1=", std::cout, EMathematicaInput);
    mat2 += mat1;
    mat2.Print("m2=", std::cout, EMathematicaInput);
}