/*		PETSc.c
 *  Only for one node.
 *      This class PETSc is created to be used a handy interface
 *  to the function calls to PETSc. For each algebric equation
 *  Ax=b, one PETSc instance is needed, since this instance has
 *  the storage of these three variables.
*/
#include "solver.h"

PETSc::PETSc()
{
	x = NULL;			/* approx solution, RHS*/
	b = NULL;
  	A = NULL;            		/* linear system matrix */

  	ksp = NULL;        		/* Krylov subspace method context */
	nullsp = NULL;
	pc = NULL;

	ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);
}

PETSc::PETSc(int ilower, int iupper, int d_nz, int o_nz)
{
	x = NULL;      			/* approx solution, RHS*/
	b = NULL;
  	A = NULL;            		/* linear system matrix */

  	ksp = NULL;          		/* Krylov subspace method context */
	nullsp = NULL;
	pc = NULL;
	Create(ilower, iupper, d_nz, o_nz);
	ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);
}

void PETSc::Create(int ilower, int iupper, int d_nz, int o_nz)
{
	Create(PETSC_COMM_WORLD, ilower, iupper, d_nz, o_nz);
}

void PETSc::Create(
	MPI_Comm Comm,
	int ilower,
	int iupper,
	int d_nz,
	int o_nz)
{
	int n	= iupper - ilower +1;

	comm 	= Comm;
	iLower	= ilower;
	iUpper 	= iupper;

        // A
	MatCreateMPIAIJ(PETSC_COMM_WORLD,n,n,PETSC_DETERMINE,PETSC_DETERMINE,
				d_nz,PETSC_NULL,o_nz,PETSC_NULL,&A);
	ierr = PetscObjectSetName((PetscObject) A, "A");
	ierr = MatSetFromOptions(A);

	// b
	ierr = VecCreate(PETSC_COMM_WORLD, &b);
	ierr = PetscObjectSetName((PetscObject) b, "b");
	ierr = VecSetSizes(b, n, PETSC_DECIDE);
	ierr = VecSetFromOptions(b);

        // x
	ierr = VecCreate(PETSC_COMM_WORLD,&x);
	ierr = PetscObjectSetName((PetscObject) x, "X");
	ierr = VecSetSizes(x, n, PETSC_DECIDE);
	ierr = VecSetFromOptions(x);
}

PETSc::~PETSc()
{
        //for petsc-3.1
	if(x!=NULL)
	{
		VecDestroy(x);
		x = NULL;
	}
	if(b!=NULL)
	{
		VecDestroy(b);
		b = NULL;
	}
	if(A!=NULL)
	{
		MatDestroy(A);
		A = NULL;
	}
	if(ksp!=NULL)
	{
		KSPDestroy(ksp);
		ksp = NULL;
	}
	if(nullsp!=NULL)
	{
		MatNullSpaceDestroy(nullsp);
		nullsp = NULL;
	}

/*
        //for petsc-3.2
        VecDestroy(&x);
        VecDestroy(&b);
        MatDestroy(&A);
        KSPDestroy(&ksp);
        MatNullSpaceDestroy(&nullsp);
*/
}

void PETSc::Reset_A()	// Reset all entries to zero ;
{
	MatZeroEntries(A);
}
void PETSc::Reset_b()  //  Reset all entries to zero ;
{
        VecZeroEntries(b);
}
void PETSc::Reset_x()
{
        VecZeroEntries(x);
}

// A
void PETSc::Set_A(PetscInt i, PetscInt j, double val)	// A[i][j]=val;
{
        ierr = MatSetValues(A,1,&i,1,&j,(const PetscScalar*)&val,INSERT_VALUES);
}

void PETSc::Add_A(PetscInt i, PetscInt j, double val)	// A[i][j]+=val;
{
        ierr = MatSetValues(A,1,&i,1,&j,(const PetscScalar*)&val,ADD_VALUES);
}

void PETSc::Get_row_of_A(PetscInt i, PetscInt *ncol, PetscInt **cols, double **row)
{
	ierr = MatGetRow(A,i,ncol,(const PetscInt**)cols,
			(const PetscScalar**)row);
	ierr = MatRestoreRow(A,i,ncol,(const PetscInt**)cols,
			(const PetscScalar**)row);
}

// x
void PETSc::Set_x(PetscInt i, double val)	// x[i]=val;
{
        ierr = VecSetValues(x,1,&i,(const PetscScalar*)&val,INSERT_VALUES);
}

void PETSc::Add_x(PetscInt i, double val)	// x[i]+=val;
{
        ierr = VecSetValues(x,1,&i,(const PetscScalar*)&val,ADD_VALUES);
}

void PETSc::Set_b(PetscInt i, double val)	// x[i]=val;
{
        ierr = VecSetValues(b,1,&i,(const PetscScalar*)&val,INSERT_VALUES);
}

void PETSc::Add_b(
	PetscInt i,
	double val)	// x[i]+=val;
{
        ierr = VecSetValues(b,1,&i,(const PetscScalar*)&val,ADD_VALUES);
}

void PETSc::Get_x(double *p)
{
	PetscScalar      *values;
	VecGetArray(x,&values);
	for(int i = 0; i < iUpper-iLower+1; i++)
		p[i] = values[i];
        VecRestoreArray(x,&values);
}

void PETSc::Get_b(double *p)
{
	PetscScalar      *values;
	VecGetArray(b,&values);
	for(int i = 0; i < iUpper-iLower+1; i++)
		p[i] = values[i];
        VecRestoreArray(b,&values);
}

void PETSc::Get_x(double *p,
	int n,
	int *global_index)
{
}

void PETSc::SetMaxIter(int val)
{
	PetscInt maxits;
	PetscReal rtol, atol, dtol;

	KSPGetTolerances(ksp, &rtol, &atol, &dtol, &maxits);
	ierr = KSPSetTolerances(ksp, rtol, atol, dtol, (PetscInt)val);
}

void PETSc::SetTol(double val)
{
	PetscInt maxits;
	PetscReal rtol, atol, dtol;

	KSPGetTolerances(ksp, &rtol, &atol, &dtol, &maxits);
	ierr = KSPSetTolerances(ksp, (PetscReal)val, atol, dtol, maxits);
}

void PETSc::SetKDim(int val)
{

}

void PETSc::GetNumIterations(PetscInt *num_iterations)
{
	KSPGetIterationNumber(ksp,num_iterations);
}

void PETSc::GetFinalRelativeResidualNorm(double *rel_resid_norm)
{
	KSPGetResidualNorm(ksp,(PetscReal*)rel_resid_norm);
}

void PETSc::Solve_GMRES(void)
{
        start_clock("Before Assemble matrix and vector");
	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

  	ierr = VecAssemblyBegin(x);
  	ierr = VecAssemblyEnd(x);

  	ierr = VecAssemblyBegin(b);
  	ierr = VecAssemblyEnd(b);
	stop_clock("After Assembly matrix and vector");


        KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN);
	KSPSetType(ksp,KSPGMRES);

        KSPSetComputeSingularValues(ksp, PETSC_TRUE);
	//KSPGetPC(ksp, &pc);
	//PCSetType(pc, PCJACOBI);

        KSPSetFromOptions(ksp);
        KSPSetUp(ksp);

	start_clock("Before KSPSolve in Solve_GMRES");
        ierr = KSPSolve(ksp,b,x);
	stop_clock("After KSPSolve in Solve_GMRES");
}

void PETSc::GetExtremeSingularValues(double *max, double *min)
{
        start_clock("Before KSPComputeExtremeSingularValues");
        KSPComputeExtremeSingularValues(ksp,(PetscReal*)max,(PetscReal*)min);
        stop_clock("After KSPComputeExtremeSingularValues");
}

void PETSc::Solve(void)
{
        start_clock("Before Assemble matrix and vector");
	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

  	ierr = VecAssemblyBegin(x);
  	ierr = VecAssemblyEnd(x);

  	ierr = VecAssemblyBegin(b);
  	ierr = VecAssemblyEnd(b);
	stop_clock("After Assembly matrix and vector");


        KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN);
        KSPSetType(ksp,KSPBCGSL);
	KSPBCGSLSetEll(ksp,2);

        KSPSetComputeSingularValues(ksp, PETSC_TRUE);
	//KSPGetPC(ksp, &pc);
	//PCSetType(pc, PCJACOBI);

        KSPSetFromOptions(ksp);
        KSPSetUp(ksp);

	start_clock("Before KSPSolve");
        ierr = KSPSolve(ksp,b,x);
	stop_clock("After KSPSolve");
}

void PETSc::Solve_withPureNeumann_GMRES(void)
{
    	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

  	ierr = VecAssemblyBegin(x);
  	ierr = VecAssemblyEnd(x);

  	ierr = VecAssemblyBegin(b);
  	ierr = VecAssemblyEnd(b);

	ierr = MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,PETSC_NULL,&nullsp);
        ierr = KSPSetNullSpace(ksp,nullsp);
	ierr = MatNullSpaceRemove(nullsp,b,PETSC_NULL);


        KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN);
	KSPSetType(ksp,KSPGMRES);

        KSPSetComputeSingularValues(ksp, PETSC_TRUE);
	//KSPGetPC(ksp, &pc);
	//PCSetType(pc, PCASM);
        KSPSetFromOptions(ksp);
        KSPSetUp(ksp);
	start_clock("Before Petsc Solve in pure neumann GMRES solver");
        ierr = KSPSolve(ksp,b,x);
	stop_clock("After Petsc Solve in pure neumann GMRES solver");
}

void PETSc::Solve_withPureNeumann(void)
{
	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
    //Print_A("HZ Matrix A");

  	ierr = VecAssemblyBegin(x);
  	ierr = VecAssemblyEnd(x);

  	ierr = VecAssemblyBegin(b);
  	ierr = VecAssemblyEnd(b);
    //Print_b("HZ Vector b");

	ierr = MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,PETSC_NULL,&nullsp);
        ierr = KSPSetNullSpace(ksp,nullsp);
        //TODO && FIXME: Orthogonal with respect to null space
	ierr = MatNullSpaceRemove(nullsp,b,PETSC_NULL);

        KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN);

	//KSPSetType(ksp,KSPMINRES);
	//KSPSetType(ksp,KSPGMRES);
	//KSPSetType(ksp,KSPBCGS);
	KSPSetType(ksp,KSPBCGSL);
	KSPBCGSLSetEll(ksp,2);

        KSPSetComputeSingularValues(ksp, PETSC_TRUE);
	//KSPGetPC(ksp, &pc);
	//PCSetType(pc, PCASM);
	//PCSetType(pc, PCMG);
	//PCMGSetLevels(pc, 3, &PETSC_COMM_WORLD);
	//PCMGSetType(pc,PC_MG_MULTIPLICATIVE);
	//PCMGSetCycleType(pc,PC_MG_CYCLE_V);
	//
        KSPSetFromOptions(ksp);
        KSPSetUp(ksp);

	start_clock("Before Petsc Solve in pure neumann solver");
        ierr = KSPSolve(ksp,b,x);
	stop_clock("After Petsc Solve in pure neumann solver");
    // Print Out Matrix A and Vector b
    //Print_A("HZ Matrix A");
    //Print_b("HZ Vector b");
    //Print_x("HZ Solution x");
}

void PETSc::Print_A(const char *filename)
{
        ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
        ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
	PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_MATLAB);
        MatView(A, PETSC_VIEWER_STDOUT_WORLD);
}

void PETSc::Print_b(const char *filename)
{
        ierr = VecAssemblyBegin(b);
        ierr = VecAssemblyEnd(b);
	PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_MATLAB);
        VecView(b, PETSC_VIEWER_STDOUT_WORLD);
}

// Print out Solution
void PETSc::Print_x(const char *filename)
{
        ierr = VecAssemblyBegin(x);
        ierr = VecAssemblyEnd(x);
	PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_MATLAB);
        VecView(x, PETSC_VIEWER_STDOUT_WORLD);
}
