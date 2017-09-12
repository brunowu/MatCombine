#ifndef _MATCOMBINE_H
#define _MATCOMBINE_H

#include "petscmat.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

typedef struct _MatrixInfo{
	int n;
	int m;
	int nnz;
} MatrixInfo;

PetscErrorCode matCombine();
PetscErrorCode read_matrix_vector(Mat * A);

#endif