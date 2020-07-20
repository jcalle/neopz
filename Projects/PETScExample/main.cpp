/*T
Examples based on PETSc tutorials
*/

#include "tutorials.h"
#include <petscksp.h>
#include <chrono>
#include <pzerror.h>

static char help[] = "Solving linear system with KSP.\n\n";

int main(int argc,char **args)
{

    int ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;

    std::cout << "\n\n#######################\nIterative solver\n";
    {
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        ierr = tridiaglinsystem();
        if (ierr!=0)
        {
            DebugStop();
        }
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        std::cout << "Time difference with PETSc (iterative) = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/10e6 << "[s]" << std::endl;
    }
    // ierr = parallellinsystem(ierr);
    std::cout << "\n\n#######################\nMUMPS solver\n";
    {
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        ierr = tridiaglinsystemMUMPS();
        if (ierr!=0)
        {
            DebugStop();
        }
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        std::cout << "Time to solve with MUMPS (direct)= " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/10e6 << "[s]" << std::endl;
    }

    ierr = PetscFinalize();

    return ierr;

}