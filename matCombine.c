#include "matCombine.h"
#define PI 3.1415926

static char help[] = "Block matrix generator";

/*
PetscErrorCode read_multiple_matrix(Mat * A, PetscPrintf *num){
	PetscErrorCode ierr;
	PetscBool flaga;
	PetscViewer fd;
	PetscInt size,sizea;

	ierr=PetscOptionsGetString(NULL,PETSC_NULL,"-mfile",filea,PETSC_MAX_PATH_LEN-1,&flaga);CHKERRQ(ierr);
	if (!flaga) {
		PetscPrintf(PETSC_COMM_WORLD,"ERRO : mfile is not properly set\n",filea);
		return 0;
	}
	PetscPrintf(PETSC_COMM_WORLD,"Loading Matrix : %s\n",filea);
	ierr=PetscViewerBinaryOpen(PETSC_COMM_WORLD,filea,FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = MatCreate(PETSC_COMM_WORLD,A);CHKERRQ(ierr);
	ierr=MatLoad(*A,fd);CHKERRQ(ierr);
	ierr=PetscViewerDestroy(&fd);CHKERRQ(ierr);
	ierr=MatGetSize(*A,&size,&sizea);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"Loaded Matrix of size : %d %d\n",size,sizea);

	return 0;
}
*/

PetscErrorCode matCombine(){

	PetscInt       	ggn, ggm,k1, k2, i,m,n,k,gcol,grow, gn, gm;
	double	      	nnz, gnnz;
	PetscScalar		rnd_value, value;
	PetscInt      	start, end, rstart, rend;
	PetscRandom   	rnd;
	PetscInt 		size;
	PetscInt 		nb;
	PetscInt        ncols;
	MatInfo     	minfo,Ainfo;
	PetscBool 		nbflg, mtflg;
	PetscErrorCode 	ierr;
	Mat 			A, M;
	PetscViewer    		output_viewer;
	char           		matrixOutputFile[PETSC_MAX_PATH_LEN];
	const PetscInt    	*cols;
	const PetscScalar 	*vals;
	PetscInt		*gcols;
    PetscInt 		idxcol;
	char			type[PETSC_MAX_PATH_LEN];
	char 			*filea[PETSC_MAX_PATH_LEN];
	PetscInt        mat_num;
	PetscBool 		mat_flg;


	MPI_Comm_size(PETSC_COMM_WORLD,&size);

	PetscPrintf(PETSC_COMM_WORLD, "The MPI world size is %d \n", size);

	PetscOptionsGetStringArray(NULL,"-mat",filea,&mat_num,&mat_flg);

    MatGetSize(M,&m, &n);
    MatGetInfo(M,MAT_GLOBAL_SUM,&minfo);
	nnz = minfo.nz_used;
    PetscPrintf(PETSC_COMM_WORLD,"Already read matrix with properties : m = %d, n = %d, nnz = %g, type = real\n",m,n,nnz);
    gm = nb * m;
    gn = nb * n;
    PetscPrintf(PETSC_COMM_WORLD,"\n]>> Start the generation of matrix : gm = %d, gn = %d, type = real\n",gm,gn);

	ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  	ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,gm,gn);CHKERRQ(ierr);
  	ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  	ierr = MatSetUp(A);CHKERRQ(ierr);	
	ierr = MatGetOwnershipRange(A, &start, &end);CHKERRQ(ierr);
	ierr = MatGetOwnershipRange(M, &rstart, &rend);CHKERRQ(ierr);

	if(strcmp(type, "matline") == 0){
		for(k = 0; k< nb; k++){
    			for(i = rstart; i < rend; i++){
    				MatGetRow(M,i,&ncols,&cols,&vals);
        			PetscMalloc1(ncols,&gcols);
				idxcol = i + k*m;
					for(k1 = 0; k1 < ncols; k1++){
						gcols[k1] = cols[k1] + k*m;
					}
				ierr = MatSetValues(A,1,&idxcol,ncols,gcols,vals,INSERT_VALUES);CHKERRQ(ierr);
    				MatRestoreRow(M,i,&ncols,&cols,&vals);
	    		}
		}
		for(k2 = start; k2 < end; k2++){
			if (k2 % 100 == 0 && k2 < gm - 350){
				gcol = k;
                        	grow = k + 350;
                        	PetscRandomGetValue(rnd,&rnd_value);
                        	value=0.01*rnd_value;
                    		ierr = MatSetValues(A,1,&grow,1,&gcol,&value,INSERT_VALUES);CHKERRQ(ierr);
			}
		}
	}

        if(strcmp(type, "matblock") == 0){
                for(k = 0; k< nb; k++){
                        for(i = rstart; i < rend; i++){
                                MatGetRow(M,i,&ncols,&cols,&vals);
                                PetscMalloc1(ncols,&gcols);
                                idxcol = i + k*m;
                                        for(k1 = 0; k1 < ncols; k1++){
                                                gcols[k1] = cols[k1] + k*m;
                                        }
                                ierr = MatSetValues(A,1,&idxcol,ncols,gcols,vals,INSERT_VALUES);CHKERRQ(ierr);
                                ierr = MatSetValues(A,1,&i,ncols,gcols,vals,INSERT_VALUES);CHKERRQ(ierr);
			        MatRestoreRow(M,i,&ncols,&cols,&vals);
                        }
                }
        }
	/*Matrix assembly*/
	PetscPrintf(PETSC_COMM_WORLD,"Assembling matrix within PETSc.\n");
	MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
//	MatView(A,PETSC_VIEWER_STDOUT_WORLD);
        MatGetSize(A, &ggm, &ggn);
        MatGetInfo(A,MAT_GLOBAL_SUM,&Ainfo);
        gnnz = Ainfo.nz_used;
	PetscPrintf(PETSC_COMM_WORLD,"Finished matrix assembly.\n");
	sprintf(matrixOutputFile,"%s_nb_%d_%dx%d_%g_nnz.gz",type, nb, ggm,ggn,gnnz);
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
	MatGenbyOther();
	PetscFinalize();

	return 0;
}