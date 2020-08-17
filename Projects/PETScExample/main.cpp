/*T
Examples based on PETSc tutorials
*/

#include "tutorials.h"
#include <petscksp.h>
#include <pzerror.h>
#include "TPZPetScMatrix.h"


int main(int argc,char **args)
{
    // Just testing the class:
    static char help[] = "Testing TPZPetscMatrix";
    int ierr = PetscInitialize(&argc,&args,(char*)0,help); if (ierr) return 0;

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
    PetscFinalize();

    // Testing the systems:
    // int ierr = tridiaglinsystem(argc, args);CHKERRQ(ierr);
    // ierr = tridiaglinsystemMUMPS(argc, args);CHKERRQ(ierr);
    // ierr = parallellinsystem(argc, args);CHKERRQ(ierr);

}