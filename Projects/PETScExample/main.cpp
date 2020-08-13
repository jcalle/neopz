/*T
Examples based on PETSc tutorials
*/

#include "tutorials.h"
#include <petscksp.h>
#include <pzerror.h>


int main(int argc,char **args)
{

    int ierr = tridiaglinsystem(argc, args);CHKERRQ(ierr);

    ierr = tridiaglinsystemMUMPS(argc, args);CHKERRQ(ierr);

    ierr = parallellinsystem(argc, args);CHKERRQ(ierr);

}