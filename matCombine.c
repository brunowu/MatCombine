#include "matCombine.h"
#define PI 3.1415926

static char help[] = "Block matrix generator";


PetscErrorCode matCombine(){

	PetscInt       	k1,i,m,n;
	double	      	gnnz;
	PetscInt      	start, end;
	PetscRandom   	rnd;
	PetscInt 		size;
	PetscInt        ncols;
	MatInfo     	Ainfo;
	PetscErrorCode 	ierr;
	Mat 			A;
	PetscViewer    		output_viewer;
	char           		matrixOutputFile[PETSC_MAX_PATH_LEN];
	const PetscInt    	*cols;
	const PetscScalar 	*vals;
	PetscInt		*gcols;
	char			m1file[PETSC_MAX_PATH_LEN];
	char			m2file[PETSC_MAX_PATH_LEN];

	PetscBool 		mat1flg,mat2flg;
	PetscViewer 	m1fd,m2fd;
	PetscInt 		mat1sizea,mat1sizeb,mat2sizea,mat2sizeb;
	Mat 			mat1, mat2;
//	PetscInt 		m, n;
	PetscInt 		ind;
 
	MPI_Comm_size(PETSC_COMM_WORLD,&size);

	PetscPrintf(PETSC_COMM_WORLD, "The MPI world size is %d \n", size);

	ierr=PetscOptionsGetString(NULL,PETSC_NULL,"-m1file",m1file,PETSC_MAX_PATH_LEN-1,&mat1flg);CHKERRQ(ierr);
	if (!mat1flg) {
		PetscPrintf(PETSC_COMM_WORLD,"ERROR : The first matrix mfile is not properly set\n",m1file);
		return 0;
	}

	PetscPrintf(PETSC_COMM_WORLD,"Loading First Matrix : %s\n",m1file);
	ierr=PetscViewerBinaryOpen(PETSC_COMM_WORLD,m1file,FILE_MODE_READ,&m1fd);CHKERRQ(ierr);
	ierr = MatCreate(PETSC_COMM_WORLD,&mat1);CHKERRQ(ierr);
	ierr=MatLoad(mat1,m1fd);CHKERRQ(ierr);
	ierr=PetscViewerDestroy(&m1fd);CHKERRQ(ierr);
	ierr=MatGetSize(mat1,&mat1sizea,&mat1sizeb);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"Loaded First Matrix of size : %d %d\n",mat1sizea,mat1sizeb);

	ierr=PetscOptionsGetString(NULL,PETSC_NULL,"-m2file",m2file,PETSC_MAX_PATH_LEN-1,&mat2flg);CHKERRQ(ierr);
	if (!mat2flg) {
		PetscPrintf(PETSC_COMM_WORLD,"ERROR : The Second matrix mfile is not properly set\n",m2file);
		return 0;
	}

	PetscPrintf(PETSC_COMM_WORLD,"Loading Second Matrix : %s\n",m2file);
	ierr=PetscViewerBinaryOpen(PETSC_COMM_WORLD,m2file,FILE_MODE_READ,&m2fd);CHKERRQ(ierr);
	ierr = MatCreate(PETSC_COMM_WORLD,&mat2);CHKERRQ(ierr);
	ierr=MatLoad(mat2,m2fd);CHKERRQ(ierr);
	ierr=PetscViewerDestroy(&m2fd);CHKERRQ(ierr);
	ierr=MatGetSize(mat2,&mat2sizea,&mat2sizeb);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"Loaded Second Matrix of size : %d %d\n",mat1sizea,mat1sizeb);

	m = mat1sizea + mat2sizea;
	n = mat1sizeb + mat2sizea;

    PetscPrintf(PETSC_COMM_WORLD,"\n]>> Start the generation of matrix : m = %d, n = %d, type = real\n",m,n);

	ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  	ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,m,n);CHKERRQ(ierr);
  	ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  	ierr = MatSetUp(A);CHKERRQ(ierr);	
	ierr = MatGetOwnershipRange(A, &start, &end);CHKERRQ(ierr);

	PetscInt m1start, m1end, m2start, m2end;
        ierr = MatGetOwnershipRange(mat1, &m1start, &m1end);CHKERRQ(ierr);
        ierr = MatGetOwnershipRange(mat2, &m2start, &m2end);CHKERRQ(ierr);

	for(i = m1start; i < m1end; i++){
		MatGetRow(mat1,i,&ncols,&cols,&vals);
		//PetscPrintf(PETSC_COMM_WORLD, "val[2] = %f\n", vals[2]);
//	        PetscMalloc1(ncols,&gcols);
                ierr = MatSetValues(A,1,&i,ncols,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
                MatRestoreRow(mat1,i,&ncols,&cols,&vals);
	}

        for(i = m2start; i < m2end; i++){
                MatGetRow(mat2,i,&ncols,&cols,&vals);
                PetscMalloc1(ncols,&gcols);
		ind = i + mat1sizea;
		for(k1 = 0; k1 < ncols; k1++){
			gcols[k1] = cols[k1] + mat1sizeb;
		}
                ierr = MatSetValues(A,1,&ind,ncols,gcols,vals,INSERT_VALUES);CHKERRQ(ierr);
                MatRestoreRow(mat2,i,&ncols,&cols,&vals);
        }

	/*Matrix assembly*/
	PetscPrintf(PETSC_COMM_WORLD,"Assembling matrix within PETSc.\n");
	MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
//	MatView(A,PETSC_VIEWER_STDOUT_WORLD);
   	MatGetInfo(A,MAT_GLOBAL_SUM,&Ainfo);
//        MatGetInfo(mat2,MAT_GLOBAL_SUM,&mat2info);
    	gnnz = Ainfo.nz_used;
//	PetscPrintf(PETSC_COMM_WORLD, "mat1info = %g, mat2info = %g \n", mat1info.nz_used,mat2info.nz_used);
	PetscPrintf(PETSC_COMM_WORLD,"Finished matrix assembly.\n");
	sprintf(matrixOutputFile,"Combine_nb_%dx%d_%g_nnz.gz", m,n,gnnz);
	PetscPrintf(PETSC_COMM_WORLD,"Dumping matrix to PETSc binary %s\n",matrixOutputFile);	
	PetscViewerBinaryOpen(PETSC_COMM_WORLD,matrixOutputFile,FILE_MODE_WRITE,&output_viewer);
	PetscViewerPushFormat(output_viewer,PETSC_VIEWER_ASCII_INFO_DETAIL);
	MatView(A,output_viewer);
	PetscViewerDestroy(&output_viewer);
	PetscRandomDestroy(&rnd);
	PetscPrintf(PETSC_COMM_WORLD,"Matrix %s Dumped\n",matrixOutputFile);
		
	return ierr;
}

int main(int argc, char ** argv){

	PetscInitialize(&argc,&argv,(char *)0,help);
	matCombine();
	PetscFinalize();

	return 0;
}
