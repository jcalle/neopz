/*T
Examples based on PETSc tutorials
*/

#include "tutorials.h"
#include <petscksp.h>
#include <pzerror.h>
#include "TPZPetScMatrix.h"


int main(int argc,char **args)
{

    static char help[] = "Testing TPZPetscMatrix";
    int ierr = PetscInitialize(&argc,&args,(char*)0,help); if (ierr) return 0;
    // int ierr = tridiaglinsystem(argc, args);CHKERRQ(ierr);

    // ierr = tridiaglinsystemMUMPS(argc, args);CHKERRQ(ierr);

    // ierr = parallellinsystem(argc, args);CHKERRQ(ierr);

    // TPZPetScMatrix<REAL> mat();
    // TPZPetScMatrix<REAL> mat2(2,2);
    TPZPetScMatrix<REAL> mat3(2,2,4);
    // TPZPetScMatrix<REAL> mat4(2,2,TPZPetScMatrix<REAL>::EDense);

    mat3.Print("m=",std::cout,EMathematicaInput);
    PetscFinalize();

}