#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <malloc.h>
#include <mkl.h>

void PrintVector(int *v, int size)
{
	int i;
	for (i=0; i<size; i++) printf("%d ", v[i]);
	printf("\n");
}

void PrintMatrix(double *A, int size)
{
	int i, j, n = (int)sqrt(size);
	for(i=0; i<n; i++)
	{
		for(j=0; j<n; j++)
		{
			printf("%f ", A[j*n+i]);
		}
		printf("\n");
	}
}

void ObtainMatrixData(char *file_name, int *n_levels, int **levels_sizes, int *n_blocks, int **block_indexes, double **A)
{
	int line, i, j;
	FILE *fp;
 	char *buffer,* buff, * element;
 	const char delimiter[2] = "\n,";
 	
	fp = fopen(file_name, "rb");
	buffer = malloc(1000*sizeof(char));
 	line = 0;

 	while (feof(fp)==0)
 	{
 		line ++; i = 0;
	    fgets(buffer, 1000, fp);
	    buff = strdup(buffer);
	    element = strsep(&buff, delimiter);

	    switch(line)
	    {
	    	case 1:
	    		*n_levels = atoi(element);
	    		*levels_sizes = NULL;
	    		*levels_sizes = (int *) malloc( *n_levels * sizeof(int));
	    		break;

    		case 2:
	    		while (element != NULL)
			    {
				    (*levels_sizes)[i] = atoi(element); 
				    i++;
				    element = strsep(&buff, delimiter);
			    }
    			break;

			case 3:
	    		*n_blocks = atoi(element);
	    		*block_indexes = NULL;
				*block_indexes = (int *) malloc( *n_blocks * sizeof(int));
				*A = NULL;
				*A = (double *) malloc( (*levels_sizes)[0] * (*levels_sizes)[0] * sizeof(double));
	    		break;

			case 4:
	    		while (element != NULL)
			    {
				    (*block_indexes)[i] = atoi(element); i++;
				    element = strsep(&buff, delimiter);
			    }
    			break;

			default:
				while (element != NULL)
			    {
			    	//printf("%s to atof: %f next buff = %s\n", element, atof(element), buff);
				    (*A)[i] = atof(element); i++;
				    element = strsep(&buff, delimiter);
			    }
    			break;
	    }
	}
 	fclose (fp);    
    free(buffer);
}

void GetBlocksIndexes( int **levels_start_i, int **levels_end_i, int n_levels, int n_blocks, int *block_indexes)
{
	int i, j;
	*levels_start_i = (int *) malloc(n_levels * sizeof(int));
    *levels_end_i = (int *) malloc(n_levels * sizeof(int));
    (*levels_start_i)[0] = 0;
    j = 1;
    for (i=1; i<n_blocks; i++)
    {
    	if(block_indexes[i] == 0)
    	{
    		(*levels_start_i)[j] = i;
    		(*levels_end_i)[j-1] = i-1;
    		j++;
    	}
    }
    (*levels_end_i)[j-1] = i-1;
}

void GetIndexes(int **indexes, int start, int size, int big_size)
{
	int i, j, arrayIndex;
	arrayIndex = 0;
	*indexes = (int *) malloc((size*size) * sizeof(int));
	for (i=0; i<size; i++)
	{  
		j = start + big_size*i;
		while (j <= start + big_size*i + size -1)
		{
			(*indexes)[arrayIndex] = j;
			j++;
			arrayIndex++;
		}
	}
}

void ApplyPivot(double **A, int *ipiv, int n, int m, int resto)
{
	int i, j;
	double aux;
	for(i=0; i<n; i++)
	{
		if((ipiv[i]-1) != i)
		{
			for(j=0; j<m; j++)
			{
				aux = (*A)[(i+resto)+m*j];
				(*A)[(i+resto)+m*j] = (*A)[(ipiv[i]-1)+m*j];
				(*A)[(ipiv[i]-1)+m*j] = aux;
			}
		}
	}
}

void ReversePivot(double **A, int *ipiv, int n, int m, int resto)
{
	int i, j;
	double aux;
	for(i=n-1; i>=0; i--)
	{
		if((ipiv[i]-1) != i)
		{
			for(j=0; j<m; j++)
			{
				aux = (*A)[(i+resto)+m*j];
				(*A)[(i+resto)+m*j] = (*A)[(ipiv[i]-1)+m*j];
				(*A)[(ipiv[i]-1)+m*j] = aux;
			}
		}
	}
}

void TRSM_row( int first_block_start_i, int size_bl_start, int bl_start_lvl, int *first_block_indexes, int n_levels, int *levels_sizes, int *block_indexes, double *A, double **LU, int *levels_start_i, int *levels_end_i, int *TRSM_performed )
{
	printf("\tTRSM(s) (row) from block starting in index %d...\n", first_block_indexes[0]);
	int i, j, k, next_in_row, next_row_start, *trsm_bl_indexes, N, NRHS, LDA, INFO;
	double *aux_A, *aux_L;

	// Lapack variables
	char UPLO = 'L', TRANS = 'N', DIAG = 'U'; 	// It is unitdiagonal
	N = size_bl_start; NRHS = N; LDA = N; INFO = 0;

	next_in_row = 1;
	aux_A = malloc(size_bl_start * size_bl_start * sizeof(double));
	aux_L = malloc(size_bl_start * size_bl_start * sizeof(double));
	
	for (i=1; i<levels_end_i[bl_start_lvl-1]-levels_start_i[bl_start_lvl-1]+1; i++) // All blocks of this size except 1st
	{
		if(i==1)
			for (j=0; j<size_bl_start; j++)
				for (k=j+1; k<size_bl_start; k++)
					aux_L[ j*size_bl_start + k ] = (*LU)[ first_block_indexes[ j*size_bl_start + k ] ];
		next_row_start = first_block_start_i + levels_sizes[0] * size_bl_start * next_in_row;
		if(block_indexes[levels_start_i[bl_start_lvl-1]+i] == next_row_start)
		{
			next_in_row ++;
			GetIndexes(&trsm_bl_indexes, next_row_start, size_bl_start, levels_sizes[0]);
			printf("\t\t> TRSM number %d, of block starting in index %d\n", next_in_row-1, trsm_bl_indexes[0]);
			//PrintVector(trsm_bl_indexes, size_bl_start*size_bl_start);
			for (j=0; j<size_bl_start; j++)
				for (k=0; k<size_bl_start; k++)
					aux_A[ j*size_bl_start + k ] = A[trsm_bl_indexes[ j*size_bl_start + k ]];
			dtrtrs_(&UPLO, &TRANS, &DIAG, &N, &NRHS, aux_L, &LDA, aux_A, &LDA, &INFO);
			for (j=0; j<size_bl_start*size_bl_start; j++) (*LU)[trsm_bl_indexes[j]] = aux_A[j];
			free(trsm_bl_indexes);
		}
	}
		
	printf("\t... end of TRSM (row)\n" );
	
	*TRSM_performed = next_in_row - 1;
	free(aux_A); free(aux_L);
}

void TRSM_col( int first_block_start_i, int size_bl_start, int bl_start_lvl, int *first_block_indexes, int n_levels, int *levels_sizes, int *block_indexes, double *A, double **LU, int *levels_start_i, int *levels_end_i, int *TRSM_performed )
{
	printf("\tTRSM(s) (col) from block starting in index %d...\n", first_block_indexes[0] );
	int i, j, k, next_in_col, next_col_start, *trsm_bl_indexes, N, NRHS, LDA, INFO;
	double *aux_A, *aux_U;

	// Lapack variables
	char UPLO = 'L', TRANS = 'N', DIAG = 'N';
	N = size_bl_start; NRHS = N; LDA = N; INFO = 0;

	next_in_col = 1;
	aux_A = malloc(size_bl_start * size_bl_start * sizeof(double));
	aux_U = malloc(size_bl_start * size_bl_start * sizeof(double));
	
	for (i=1; i<levels_end_i[bl_start_lvl-1]-levels_start_i[bl_start_lvl-1]+1; i++) // All blocks of this size except 1st
	{
		if(i==1)
			for (j=0; j<size_bl_start; j++)
				for (k=j; k<size_bl_start; k++)
					aux_U[ j * size_bl_start + k ] = (*LU)[first_block_indexes[ k * size_bl_start + j ]];
		next_col_start = first_block_start_i + size_bl_start * next_in_col;
		if(block_indexes[levels_start_i[bl_start_lvl-1]+i] == next_col_start)
		{
			next_in_col ++;
			GetIndexes(&trsm_bl_indexes, next_col_start, size_bl_start, levels_sizes[0]);
			printf("\t\t> TRSM number %d, of block starting in index %d\n", next_in_col-1, trsm_bl_indexes[0]);
			//PrintVector(trsm_bl_indexes, size_bl_start*size_bl_start);
			for (j=0; j<size_bl_start; j++)
				for (k=0; k<size_bl_start; k++) // Store U matrix fragment (transposed)
					aux_A[ j * size_bl_start + k ] = A[trsm_bl_indexes[ k * size_bl_start + j ]];
			//PrintVector(first_block_indexes, size_bl_start*size_bl_start);
			dtrtrs_(&UPLO, &TRANS, &DIAG, &N, &NRHS, aux_U, &LDA, aux_A, &LDA, &INFO);
			for (k=0; k<size_bl_start; k++) // Store res in LU (transpose way)
				for (j=0; j<size_bl_start; j++)
					(*LU)[trsm_bl_indexes[ k * size_bl_start + j ]] = aux_A[ j * size_bl_start + k ];
			free(trsm_bl_indexes);
		}
	}
	printf("\t... end of TRSM (col)\n" );
	*TRSM_performed = next_in_col - 1;
	free(aux_A); free(aux_U);
}

void Update( int first_block_start_i, int size_bl_start, int bl_start_lvl, int n_levels, int *levels_sizes, int *block_indexes, double **A, double *LU, int *levels_start_i, int *levels_end_i, int to_upd )
{
	printf("\tUPD from block starting in index %d...\n", first_block_start_i);
	int bl, *upd_bl_indexes, diag_start, diag_start_row, diag_start_col, i, j, *U_idx, *L_idx, N, LDA;
	double *aux_U, *aux_L, *aux_X;

	// VARS for UPD
	N = size_bl_start;
	LDA = N;
	char TRANS = 'N';
	double ALPHA = 1.0, BETA = 0.0;

	aux_L = malloc(size_bl_start * size_bl_start * sizeof(double));
	aux_U = malloc(size_bl_start * size_bl_start * sizeof(double));
	aux_X = malloc(size_bl_start * size_bl_start * sizeof(double));

	// UPD from the first smaller block to same-sized blocks in the same column
	if(to_upd > 0)
	{
		
		for(bl=1; bl<=to_upd; bl++)
		{
			i = first_block_start_i + levels_sizes[0] * size_bl_start * bl;
			GetIndexes(&U_idx, i, size_bl_start, levels_sizes[0]);
			i = first_block_start_i + size_bl_start * bl;
			GetIndexes(&L_idx, i, size_bl_start, levels_sizes[0]);
			for (j=0; j<size_bl_start*size_bl_start; j++)
			{
				aux_U[j] = LU[ U_idx[j] ];
				aux_L[j] = LU[ L_idx[j] ];
			}
			diag_start = first_block_start_i + bl * levels_sizes[0] * size_bl_start + bl * size_bl_start;
			GetIndexes(&upd_bl_indexes, diag_start, size_bl_start, levels_sizes[0]);
			printf("\t\t> Updating block starting in index %d\n", upd_bl_indexes[0]);
			//PrintVector(upd_bl_indexes, size_bl_start*size_bl_start);
			dgemm_(&TRANS, &TRANS, &N, &N, &N, &ALPHA, aux_L, &LDA, aux_U, &LDA, &BETA, aux_X, &N);
			for (i=0; i<size_bl_start*size_bl_start; i++) (*A)[upd_bl_indexes[i]] = (*A)[upd_bl_indexes[i]] - aux_X[i];
			for(j=1; j<=to_upd-bl; j++) // Blocks on the right side of diagonal block (in the row)
			{
				printf("diag_start = %d\n", diag_start);
				printf("U_idx = "); PrintVector(U_idx, size_bl_start*size_bl_start);
				printf("L_idx = "); PrintVector(L_idx, size_bl_start*size_bl_start);
				i = first_block_start_i + levels_sizes[0] * size_bl_start * (j+1);
				GetIndexes(&U_idx, i, size_bl_start, levels_sizes[0]);
				for (i=0; i<size_bl_start*size_bl_start; i++) aux_U[i] = LU[ U_idx[i] ];
				diag_start_row = diag_start + size_bl_start * levels_sizes[0] * j;
				GetIndexes(&upd_bl_indexes, diag_start_row, size_bl_start, levels_sizes[0]);
				printf("\t\t> Updating block starting in index %d\n", diag_start_row);
				//PrintVector(upd_bl_indexes, size_bl_start*size_bl_start);
				dgemm_(&TRANS, &TRANS, &N, &N, &N, &ALPHA, aux_L, &LDA, aux_U, &LDA, &BETA, aux_X, &N);
				for (i=0; i<size_bl_start*size_bl_start; i++) (*A)[upd_bl_indexes[i]] = (*A)[upd_bl_indexes[i]] - aux_X[i];
			}
			i = first_block_start_i + levels_sizes[0] * size_bl_start * bl;
			GetIndexes(&U_idx, i, size_bl_start, levels_sizes[0]);
			for (j=0; j<size_bl_start*size_bl_start; j++) aux_U[j] = LU[ U_idx[j] ];
			for(j=1; j<=to_upd-bl; j++) // Blocks under diagonal block (in the column)
			{
				printf("diag_start = %d\n", diag_start);
				printf("U_idx = "); PrintVector(U_idx, size_bl_start*size_bl_start);
				printf("L_idx = "); PrintVector(L_idx, size_bl_start*size_bl_start);
				i = first_block_start_i + size_bl_start * (j+1);
				GetIndexes(&L_idx, i, size_bl_start, levels_sizes[0]);
				for (i=0; i<size_bl_start*size_bl_start; i++) aux_L[i] = LU[ L_idx[i] ];
				diag_start_col = diag_start + size_bl_start * j;
				GetIndexes(&upd_bl_indexes, diag_start_col, size_bl_start, levels_sizes[0]);
				printf("\t\t> Updating block starting in index %d\n", diag_start_col);
				//PrintVector(upd_bl_indexes, size_bl_start*size_bl_start);
				dgemm_(&TRANS, &TRANS, &N, &N, &N, &ALPHA, aux_L, &LDA, aux_U, &LDA, &BETA, aux_X, &N);
				for (i=0; i<size_bl_start*size_bl_start; i++) (*A)[upd_bl_indexes[i]] = (*A)[upd_bl_indexes[i]] - aux_X[i];
			}
		}
	 	free(aux_X); free(U_idx); free(L_idx); free(upd_bl_indexes);
	}
	printf("\t... end of UPD\n");
	free(aux_U); free(aux_L);
}

void CheckLU(double *LU, int size, int *BIG_IPIV)
{
	int i,j;
	char SIDE = 'L', UPLO = 'L', TRANSA = 'N', DIAG = 'U';
	double ALPHA = 1.0, *U = malloc(size*size*sizeof(double));
	for(i=0; i<size; i++)
		for(j=0;j<size;j++)
		{
			if(j>i) U[i*size+j] = 0.0;
			else U[i*size+j] = LU[i*size+j];
		}
	dtrmm_(&SIDE, &UPLO, &TRANSA, &DIAG, &size, &size, &ALPHA, LU, &size, U, &size);
	ReversePivot(&U, BIG_IPIV, size, size, 0);
	PrintMatrix(U, size*size);
	free(U);
}

int main (int argc, char *argv[])
{

	// ----------------------- CHECK USAGE 
	if (argc != 2)
	{
	    printf ("Usage: Iterative_LU <FileOfMatrix>\n"); return 0;
	}

	// ----------------------- VARIABLES DECLATIONS 
	int i, j, k, n_levels, n_blocks, A_size, *block_indexes, *levels_sizes, *levels_start_i, *levels_end_i;
	int first_block_start_i, *first_block_indexes, n_blocks_level, N, LDA, INFO, NRHS;
	int next_in_row, next_in_col, next_row_start, next_col_start;
	int diag_start, *U_idx, *L_idx, to_upd, i_diag, diag_start_i, *diag_bl_indexes;
	int prev, keep_prev, previous_lvl_bl_start_i, *previous_lvl_bl_indexes, TRSMs_in_row, TRSMs_in_col;
	double *A, *LU, *aux_A;
	char *file_name;
	
	// ----------------------- GET DATA FROM INPUT
	file_name = argv[1];

	// ----------------------- READ MATRIX DATA FROM FILE 
	printf("\nGiven file containing matrix information: %s\n\n", file_name);
 	printf("Obtaining matrix information... ");
 	ObtainMatrixData(file_name, &n_levels, &levels_sizes, &n_blocks, &block_indexes, &A);
 	printf("DONE\n\n");
 	
	// ----------------------- GET LEVELS INDEXES
    printf("Obtaining levels indexes... ");
    GetBlocksIndexes(&levels_start_i, &levels_end_i, n_levels, n_blocks, block_indexes); 
    printf("DONE\n\n");
    
	// ----------------------- PERFORM LU    
    printf("Computing L and U matrices... \n");

    // Initialize L, U matrices to zeros
	LU = malloc(levels_sizes[0] * levels_sizes[0] * sizeof(double));
	aux_A = malloc(levels_sizes[n_levels-1]*levels_sizes[n_levels-1] * sizeof(double));
				
	// Vars for LU 
	INFO = 0;
    N = levels_sizes[n_levels-1];
    LDA = levels_sizes[n_levels-1];
    int IPIV[LDA]; int BIG_IPIV[levels_sizes[0]];

    // LU of the first smallest block
	first_block_start_i = 0;
	n_blocks_level = 1;
	first_block_indexes = malloc(n_blocks_level * sizeof(int));
	GetIndexes(&first_block_indexes, first_block_start_i, levels_sizes[n_levels-1], levels_sizes[0]);
	printf("\n\tLU of block starting in index %d...\n", first_block_indexes[0]);
	for (i=0; i<levels_sizes[n_levels-1]*levels_sizes[n_levels-1]; i++) aux_A[i] = A[first_block_indexes[i]];
	dgetrf_(&N, &N, aux_A, &LDA, IPIV, &INFO);
	//PrintVector(IPIV, LDA);
	ApplyPivot(&A, IPIV, LDA, levels_sizes[0], 0);
	for (i=0; i< levels_sizes[0]; i++) BIG_IPIV[i] = IPIV[i];
	for (i=0; i<levels_sizes[n_levels-1]*levels_sizes[n_levels-1]; i++) LU[first_block_indexes[i]] = aux_A[i]; 
	printf("\t... end of LU\n" );

	// Vars for TRSM (row)
	INFO = 0;
    N = levels_sizes[n_levels-1];
    LDA = N;
    NRHS = N;

    // TRSM from the first smaller block to same-sized blocks in the same row
	TRSM_row( first_block_start_i, levels_sizes[n_levels-1], n_levels, first_block_indexes, n_levels, levels_sizes, block_indexes, A, &LU, levels_start_i, levels_end_i, &TRSMs_in_row );
	
	// TRSM from the first smaller block to same-sized blocks in the same column
	TRSM_col( first_block_start_i, levels_sizes[n_levels-1], n_levels, first_block_indexes, n_levels, levels_sizes, block_indexes, A, &LU, levels_start_i, levels_end_i, &TRSMs_in_col );

	// UPD from the first smaller block to same-sized blocks in the diagonal
	if ((TRSMs_in_row) < (TRSMs_in_col)) to_upd = TRSMs_in_row;
	else to_upd = TRSMs_in_col;
	Update( first_block_start_i, levels_sizes[n_levels-1], n_levels, n_levels, levels_sizes, block_indexes, &A, LU, levels_start_i, levels_end_i, to_upd );

	// REMARK: diagonal is made of smallest level blocks
	// Move along diagonal blocks and execute: LU + TRSM(s) + UPD(s)
	for(i_diag=1; i_diag<levels_sizes[0]/levels_sizes[n_levels-1]; i_diag++)
	{
		diag_start_i = i_diag * levels_sizes[0] * levels_sizes[n_levels-1] + i_diag * levels_sizes[n_levels-1];
		GetIndexes(&diag_bl_indexes, diag_start_i, levels_sizes[n_levels-1], levels_sizes[0]);
		
		// LU of diagonal block
		printf("\n\tLU of block starting in index %d...\n", diag_start_i);
		for (i=0; i<levels_sizes[n_levels-1]*levels_sizes[n_levels-1]; i++) aux_A[i] = A[diag_bl_indexes[i]];
		dgetrf_(&N, &N, aux_A, &LDA, IPIV, &INFO);
		j = diag_bl_indexes[0] % levels_sizes[0];
		for (i=0; i<LDA; i++) { IPIV[i] = IPIV[i] + j; BIG_IPIV[i + j] = IPIV[i]; }
		//PrintVector(IPIV, LDA);
		ApplyPivot(&A, IPIV, LDA, levels_sizes[0], j);
		ApplyPivot(&LU, IPIV, LDA, levels_sizes[0], j);
		for (i=0; i<levels_sizes[n_levels-1]*levels_sizes[n_levels-1]; i++) LU[diag_bl_indexes[i]] = aux_A[i]; 
		printf("\t... end of LU\n");

		// TRSM of next blocks of same level in same row that diagonal block
		TRSM_row( diag_start_i, levels_sizes[n_levels-1], n_levels, diag_bl_indexes, n_levels, levels_sizes, block_indexes, A, &LU, levels_start_i, levels_end_i, &TRSMs_in_row );

		// TRSM of next blocks of same level in same column that diagonal block
		TRSM_col( diag_start_i, levels_sizes[n_levels-1], n_levels, diag_bl_indexes, n_levels, levels_sizes, block_indexes, A, &LU, levels_start_i, levels_end_i, &TRSMs_in_col );
		
		// UPD of next diagonal blocks of same level (as much as TRSMs have been performed)
		// If no TRSM's have been performed, then evaluate if there is a bigger block to perform TRSM(s) and UPD(s)
		if ((TRSMs_in_row) < (TRSMs_in_col)) to_upd = TRSMs_in_row;
		else to_upd = TRSMs_in_col;
		if(to_upd > 0)
		{
			Update( diag_start_i, levels_sizes[n_levels-1], n_levels-1, n_levels, levels_sizes, block_indexes, &A, LU, levels_start_i, levels_end_i, to_upd );
		}
		else
		{
			prev = 1;
			keep_prev = 1;
			previous_lvl_bl_start_i = 2;
			while (keep_prev && (previous_lvl_bl_start_i > 0))
			{
				previous_lvl_bl_start_i = diag_bl_indexes[levels_sizes[n_levels-1]*levels_sizes[n_levels-1]-1] - ( (levels_sizes[n_levels-1-prev] - 1) *levels_sizes[0] + levels_sizes[n_levels-1-prev] - 1);
				GetIndexes(&previous_lvl_bl_indexes, previous_lvl_bl_start_i, levels_sizes[n_levels-1-prev], levels_sizes[0]);
				printf("\tNo TRSM(s) have been performed --> Try TRSM(s) of bigger levels blocks.\n\tPrev level block start i: %d (size = %dx%d)\n", previous_lvl_bl_start_i, levels_sizes[n_levels-1-prev], levels_sizes[n_levels-1-prev]);
				
				// TRSM of next blocks of same level in same row that bigger block
				TRSM_row( previous_lvl_bl_start_i, levels_sizes[n_levels-1-prev], n_levels-prev, previous_lvl_bl_indexes, n_levels, levels_sizes, block_indexes, A, &LU, levels_start_i, levels_end_i, &TRSMs_in_row );

				// TRSM of next blocks of same level in same column that bigger block
				TRSM_col( previous_lvl_bl_start_i, levels_sizes[n_levels-1-prev], n_levels-prev, previous_lvl_bl_indexes, n_levels, levels_sizes, block_indexes, A, &LU, levels_start_i, levels_end_i, &TRSMs_in_col );

				// UPD of next diagonal blocks of same level (as much as TRSMs have been performed)
				if ((TRSMs_in_row) < (TRSMs_in_col)) to_upd = TRSMs_in_row;
				else to_upd = TRSMs_in_col;
				if(to_upd > 0)
				{
					keep_prev = 0;
					Update( previous_lvl_bl_start_i, levels_sizes[n_levels-1-prev], n_levels-prev-1, n_levels, levels_sizes, block_indexes, &A, LU, levels_start_i, levels_end_i, to_upd );
				}
				else prev++;
				free(previous_lvl_bl_indexes);
			}
		}
		free(diag_bl_indexes);
	}
	
	printf("\nDONE\n\n");

	// Print solution (L and U matrices have been stored in the same matrix)
	printf("Computed matrix LU is: \n\n");
	PrintMatrix(LU, levels_sizes[0] * levels_sizes[0]);
	printf("\nL * U =\n\n");
 	CheckLU(LU, levels_sizes[0], BIG_IPIV);

    // ----------------------- FREE MEMORY
    printf("\nFree memory... ");
    free(A); free(LU); free(aux_A);
   	free(levels_sizes); free(block_indexes);
    free(levels_start_i); free(levels_end_i);
    free(first_block_indexes);
    printf("DONE\n\n");
}
