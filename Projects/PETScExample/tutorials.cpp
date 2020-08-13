#include "tutorials.h"
#include <petscksp.h>
#include <chrono>
#include "pzerror.h"

void testingpetsc(int argc,char **args)
{
    static char help[] = "Solving linear system with KSP.\n\n";
    int ierr = PetscInitialize(&argc,&args,(char*)0,help); if (ierr) return;

    std::cout << "\n\n#######################\nIterative solver\n";
    {
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        ierr = tridiaglinsystem(argc, args);
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
        ierr = tridiaglinsystemMUMPS(argc, args);
        if (ierr!=0)
        {
            DebugStop();
        }
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        std::cout << "Time to solve with MUMPS (direct)= " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/10e6 << "[s]" << std::endl;
    }

    ierr = PetscFinalize();
    
}

int tridiaglinsystem(int argc,char **args)
{
    Vec x, b, u; // approx solution, RHS, exact solution 
    Mat A; // linear system matrix */
    KSP ksp; // linear solver context */
    PC pc; // preconditioner context */
    PetscReal norm; // norm of solution error */
    PetscInt i,n = 1000000,col[3],its;
    PetscMPIInt size;
    PetscScalar value[3];
    int ierr;
    static char help[] = "Solving linear system with KSP.\n\n";

    // MPI function: Determines the size of the group associated with a communicator
    // PETSC_COMM_WORLD: communicator
    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);
    CHKERRQ(ierr);

    // Checks if the number of MPI ranks is different from 1 (this example is serial!)
    if (size != 1)
    {
        SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_WRONG_MPI_SIZE,"This is a uniprocessor example only!");
    }
    // Gets the integer value for a particular option in the database 
    // n - the integer value to return.
    // ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Compute the matrix and right-hand-side vector that define
    the linear system, Ax = b.
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /*
    Create vectors. Note that we form 1 vector from scratch and
    then duplicate as needed.
    */
    ierr = VecCreate(PETSC_COMM_WORLD,&x); CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject) x, "Solution"); CHKERRQ(ierr);
    int dim = 1000000;
    ierr = VecSetSizes(x,dim,n);CHKERRQ(ierr);
    ierr = VecSetFromOptions(x);CHKERRQ(ierr);

    ierr = VecDuplicate(x,&b);CHKERRQ(ierr);
    ierr = VecDuplicate(x,&u);CHKERRQ(ierr);

    /*
    Create matrix. When using MatCreate(), the matrix format can
    be specified at runtime.
    Performance tuning note: For problems of substantial size,
    preallocation of matrix memory is crucial for attaining good
    performance. See the matrix chapter of the users manual for details.
    */
    ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
    ierr = MatSetSizes(A,dim,dim,n,n);CHKERRQ(ierr);
    ierr = MatSetFromOptions(A);CHKERRQ(ierr);
    ierr = MatSetUp(A);CHKERRQ(ierr);

    /*
    Assemble matrix
    */
    value[0] = -1.0; value[1] = 2.0; value[2] = -1.0;
    for (i=1; i<n-1; i++) {
        col[0] = i-1; col[1] = i; col[2] = i+1;
        ierr = MatSetValues(A,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);
    }
    i = n - 1; col[0] = n - 2; col[1] = n - 1;
    ierr = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
    i = 0; col[0] = 0; col[1] = 1; value[0] = 2.0; value[1] = -1.0;
    ierr = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    /*
    Set exact solution; then compute right-hand-side vector.
    */
    ierr = VecSet(u,1.0);CHKERRQ(ierr);
    ierr = MatMult(A,u,b);CHKERRQ(ierr);
    //Print the matrix
    bool print = false;
    if(print){
        MatView(A, PETSC_VIEWER_STDOUT_SELF);
    }
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Create the linear solver and set various options
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
    /*
    Set operators. Here the matrix that defines the linear system
    also serves as the matrix that defines the preconditioner.
    */
    ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
    /*
    Set linear solver defaults for this problem (optional).
    - By extracting the KSP and PC contexts from the KSP context,
    we can then directly call any KSP and PC routines to set
    various options.
    - The following four statements are optional; all of these
    parameters could alternatively be specified at runtime via
    KSPSetFromOptions();
    */
    ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
    ierr = PCSetType(pc,PCLU);CHKERRQ(ierr);
    ierr = KSPSetTolerances(ksp,1.e-8,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
    CHKERRQ(ierr);
    /*
    Set runtime options, e.g.,
    -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
    These options will override those specified above as long as
    KSPSetFromOptions() is called _after_ any other customization
    routines.
    */

    ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Solve the linear system
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
    /*
    View solver info; we could instead use the option -ksp_view to
    print this info to the screen at the conclusion of KSPSolve().
    */
    ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Check the solution and clean up
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    ierr = VecAXPY(x,-1.0,u);CHKERRQ(ierr);
    ierr = VecNorm(x,NORM_2,&norm);CHKERRQ(ierr);
    ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of error %g, Iterations %D\n",(double)norm,its);CHKERRQ(ierr);
    /*
    Free work space. All PETSc objects should be destroyed when they
    are no longer needed.
    */
    ierr = VecDestroy(&x);CHKERRQ(ierr); ierr = VecDestroy(&u);CHKERRQ(ierr);
    ierr = VecDestroy(&b);CHKERRQ(ierr); ierr = MatDestroy(&A);CHKERRQ(ierr);
    ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
    
    return ierr;
}

int tridiaglinsystemMUMPS(int argc,char **args)
{
    static char help[] = "Solving linear system with MUMPS.\n\n";
    int ierr;
    Vec x, b, u; // approx solution, RHS, exact solution 
    Mat A, F; // linear system matrix */
    KSP ksp; // linear solver context */
    PC pc; // preconditioner context */
    PetscReal norm; // norm of solution error */
    PetscInt i,n = 1000000,col[3],its;
    PetscMPIInt size;
    PetscScalar value[3];

    // MPI function: Determines the size of the group associated with a communicator
    // PETSC_COMM_WORLD: communicator
    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);
    CHKERRQ(ierr);

    // Checks if the number of MPI ranks is different from 1 (this example is serial!)
    if (size != 1)
    {
        SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_WRONG_MPI_SIZE,"This is a uniprocessor example only!");
    }
    // Gets the integer value for a particular option in the database 
    // n - the integer value to return.
    // ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Compute the matrix and right-hand-side vector that define
    the linear system, Ax = b.
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /*
    Create vectors. Note that we form 1 vector from scratch and
    then duplicate as needed.
    */
    ierr = VecCreate(PETSC_COMM_WORLD,&x); CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject) x, "Solution"); CHKERRQ(ierr);
    int dim = 1000000;
    ierr = VecSetSizes(x,dim,n);CHKERRQ(ierr);
    ierr = VecSetFromOptions(x);CHKERRQ(ierr);

    ierr = VecDuplicate(x,&b);CHKERRQ(ierr);
    ierr = VecDuplicate(x,&u);CHKERRQ(ierr);

    /*
    Create matrix. When using MatCreate(), the matrix format can
    be specified at runtime.
    Performance tuning note: For problems of substantial size,
    preallocation of matrix memory is crucial for attaining good
    performance. See the matrix chapter of the users manual for details.
    */
    ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
    ierr = MatSetSizes(A,dim,dim,n,n);CHKERRQ(ierr);
    ierr = MatSetFromOptions(A);CHKERRQ(ierr);
    ierr = MatSetUp(A);CHKERRQ(ierr);

    /*
    Assemble matrix
    */
    value[0] = -1.0; value[1] = 2.0; value[2] = -1.0;
    for (i=1; i<n-1; i++) {
        col[0] = i-1; col[1] = i; col[2] = i+1;
        ierr = MatSetValues(A,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);
    }
    i = n - 1; col[0] = n - 2; col[1] = n - 1;
    ierr = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
    i = 0; col[0] = 0; col[1] = 1; value[0] = 2.0; value[1] = -1.0;
    ierr = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    
    //print matrix 
    bool print = 0;
    if(print){
        MatView(A, PETSC_VIEWER_STDOUT_SELF);
    }
    /*
    Set exact solution; then compute right-hand-side vector.
    */
    ierr = VecSet(u,1.0);CHKERRQ(ierr);
    ierr = MatMult(A,u,b);CHKERRQ(ierr);
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Create the linear solver and set various options
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
    /*
    Set operators. Here the matrix that defines the linear system
    also serves as the matrix that defines the preconditioner.
    */
    ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
    /*
    Set linear solver defaults for this problem (optional).
    - By extracting the KSP and PC contexts from the KSP context,
    we can then directly call any KSP and PC routines to set
    various options.
    - The following four statements are optional; all of these
    parameters could alternatively be specified at runtime via
    KSPSetFromOptions();
    */
    KSPSetOperators(ksp,A,A);
    ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
    ierr = PCSetType(pc,PCLU);CHKERRQ(ierr);
    PCFactorSetMatSolverType(pc,MATSOLVERMUMPS);
    PCFactorSetUpMatSolverType(pc);
    PCFactorGetMatrix(pc,&A);
    int icntl=7; int ival = 2;
    MatMumpsSetIcntl(A,icntl,ival);

    /*
    Set runtime options, e.g.,
    -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
    These options will override those specified above as long as
    KSPSetFromOptions() is called _after_ any other customization
    routines.
    */

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Solve the linear system
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
    /*
    View solver info; we could instead use the option -ksp_view to
    print this info to the screen at the conclusion of KSPSolve().
    */
    ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Check the solution and clean up
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    ierr = VecAXPY(x,-1.0,u);CHKERRQ(ierr);
    ierr = VecNorm(x,NORM_2,&norm);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of error %g\n",(double)norm);CHKERRQ(ierr);
    /*
    Free work space. All PETSc objects should be destroyed when they
    are no longer needed.
    */
    // ierr = VecDestroy(&x);CHKERRQ(ierr); ierr = VecDestroy(&u);CHKERRQ(ierr);
    // ierr = VecDestroy(&b);CHKERRQ(ierr); ierr = MatDestroy(&A);CHKERRQ(ierr);
    // ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
    
    return ierr;
}

int parallellinsystem(int argc,char **args)
{
    int ierr;
    Vec x,b,u; /* approx solution, RHS, exact solution */
    Mat A; /* linear system matrix */
    KSP ksp; /* linear solver context */
    PetscReal norm; /* norm of solution error */
    PetscInt i,j,Ii,J,Istart,Iend,m = 8,n = 7,its;
    PetscBool flg;
    PetscScalar v;

    ierr = PetscOptionsGetInt(NULL,NULL,"-m",&m,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Compute the matrix and right-hand-side vector that define
    the linear system, Ax = b.
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /*
    Create parallel matrix, specifying only its global dimensions.
    When using MatCreate(), the matrix format can be specified at
    runtime. Also, the parallel partitioning of the matrix is
    determined by PETSc at runtime.
    Performance tuning note: For problems of substantial size,
    preallocation of matrix memory is crucial for attaining good
    performance. See the matrix chapter of the users manual for details.
    */
    n=100;m=100;
    ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
    ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,m*n,m*n);CHKERRQ(ierr);
    ierr = MatSetFromOptions(A);CHKERRQ(ierr);
    ierr = MatMPIAIJSetPreallocation(A,5,NULL,5,NULL);CHKERRQ(ierr);
    ierr = MatSeqAIJSetPreallocation(A,5,NULL);CHKERRQ(ierr);
    ierr = MatSeqSBAIJSetPreallocation(A,1,5,NULL);CHKERRQ(ierr);
    ierr = MatMPISBAIJSetPreallocation(A,1,5,NULL,5,NULL);CHKERRQ(ierr);
    ierr = MatMPISELLSetPreallocation(A,5,NULL,5,NULL);CHKERRQ(ierr);
    ierr = MatSeqSELLSetPreallocation(A,5,NULL);CHKERRQ(ierr);
    /*
    Currently, all PETSc parallel matrix formats are partitioned by
    contiguous chunks of rows across the processors. Determine which
    rows of the matrix are locally owned.
    */
    ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
    /*
    Set matrix elements for the 2-D, five-point stencil in parallel.
    - Each processor needs to insert only elements that it owns
    locally (but any non-local elements will be sent to the
    appropriate processor during matrix assembly).
    - Always specify global rows and columns of matrix entries.
    Note: this uses the less common natural ordering that orders first
    all the unknowns for x = h then for x = 2h etc; Hence you see J = Ii +- n
    instead of J = I +- m as you might expect. The more standard ordering
    would first do all variables for y = h, then y = 2h etc.
    */

    for (Ii=Istart; Ii<Iend; Ii++) {
        v = -1.0; i = Ii/n; j = Ii - i*n;
        if (i>0) {
            J = Ii - n; ierr = MatSetValues(A,1,&Ii,1,&J,&v,ADD_VALUES);CHKERRQ(ierr);}
        if (i<m-1) {
            J = Ii + n; ierr = MatSetValues(A,1,&Ii,1,&J,&v,ADD_VALUES);CHKERRQ(ierr);}
        if (j>0) {
            J = Ii - 1; ierr = MatSetValues(A,1,&Ii,1,&J,&v,ADD_VALUES);CHKERRQ(ierr);}
        if (j<n-1) {
            J = Ii + 1; ierr = MatSetValues(A,1,&Ii,1,&J,&v,ADD_VALUES);CHKERRQ(ierr);
        }
        v = 4.0;
        ierr = MatSetValues(A,1,&Ii,1,&Ii,&v,ADD_VALUES);CHKERRQ(ierr);
    }
    /*
    Assemble matrix, using the 2-step process:
    MatAssemblyBegin(), MatAssemblyEnd()
    Computations can be done while messages are in transition
    by placing code between these two statements.
    */
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    /* A is symmetric. Set symmetric flag to enable ICC/Cholesky preconditioner */
    ierr = MatSetOption(A,MAT_SYMMETRIC,PETSC_TRUE);CHKERRQ(ierr);
    /*
    Create parallel vectors.
    - We form 1 vector from scratch and then duplicate as needed.
    - When using VecCreate(), VecSetSizes and VecSetFromOptions()
    in this example, we specify only the
    vector’s global dimension; the parallel partitioning is determined
    at runtime.
    - When solving a linear system, the vectors and matrices MUST
    be partitioned accordingly. PETSc automatically generates
    appropriately partitioned matrices and vectors when MatCreate()
    and VecCreate() are used with the same communicator.
    - The user can alternatively specify the local vector and matrix
    dimensions when more sophisticated partitioning is needed
    (replacing the PETSC DECIDE argument in the VecSetSizes() statement
    below).
    */
    ierr = VecCreate(PETSC_COMM_WORLD,&u);CHKERRQ(ierr);
    ierr = VecSetSizes(u,PETSC_DECIDE,m*n);CHKERRQ(ierr);
    ierr = VecSetFromOptions(u);CHKERRQ(ierr);
    ierr = VecDuplicate(u,&b);CHKERRQ(ierr);
    ierr = VecDuplicate(b,&x);CHKERRQ(ierr);

    /*
    Set exact solution; then compute right-hand-side vector.
    By default we use an exact solution of a vector with all
    elements of 1.0;
    */
    ierr = VecSet(u,1.0);CHKERRQ(ierr);
    ierr = MatMult(A,u,b);CHKERRQ(ierr);
    /*
    View the exact solution vector if desired
    */
    flg = PETSC_FALSE;
    ierr = PetscOptionsGetBool(NULL,NULL,"-view_exact_sol",&flg,NULL);CHKERRQ(ierr);
    if (flg) {
        ierr = VecView(u,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    }
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Create the linear solver and set various options
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
    /*
    Set operators. Here the matrix that defines the linear system
    also serves as the preconditioning matrix.
    */
    ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);

    /*
    Set linear solver defaults for this problem (optional).
    - By extracting the KSP and PC contexts from the KSP context,
    we can then directly call any KSP and PC routines to set
    various options.
    - The following two statements are optional; all of these
    parameters could alternatively be specified at runtime via
    KSPSetFromOptions(). All of these defaults can be
    overridden at runtime, as indicated below.
    */
    ierr = KSPSetTolerances(ksp,1.e-2/((m+1)*(n+1)),1.e-50,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
    /*
    Set runtime options, e.g.,
    -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
    These options will override those specified above as long as
    KSPSetFromOptions() is called _after_ any other customization
    routines.
    */
    ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Solve the linear system
    ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Check the solution and clean up
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    ierr = VecAXPY(x,-1.0,u);CHKERRQ(ierr);
    ierr = VecNorm(x,NORM_2,&norm);CHKERRQ(ierr);
    ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
    /*
    Print convergence information. PetscPrintf() produces a single
    print statement from all processes that share a communicator.
    An alternative is PetscFPrintf(), which prints to a file.
    */
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of error %g iterations %D\n",(double)norm,its);CHKERRQ(ierr);
    /*
    Free work space. All PETSc objects should be destroyed when they
    are no longer needed.
    */
    ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
    ierr = VecDestroy(&u);CHKERRQ(ierr); ierr = VecDestroy(&x);CHKERRQ(ierr);
    ierr = VecDestroy(&b);CHKERRQ(ierr); ierr = MatDestroy(&A);CHKERRQ(ierr);

    return ierr;

}