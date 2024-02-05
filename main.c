#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "functions.h"

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        fprintf(stderr, "Usage: %s <matrix_file.mtx>\n", argv[0]);
        return EXIT_FAILURE;
    }

    // <Handle the inputs here>
    const char *filename = argv[1];

    CSRMatrix A;
    ReadMMtoCSR(filename, &A);

    // Initializing the vector b and x(in Ax=b)
    double *x = (double *)malloc(A.num_cols * sizeof(double));
    double *b = (double *)malloc(A.num_rows * sizeof(double));

    // Set all elements of b to 1 and x to 0
    for (int i = 0; i < A.num_rows; ++i)
    {
        b[i] = 1.0;
        x[i] = 2.0;
    }

    //calling on the solver and recording the time it takes
    clock_t start_time = clock();
    solver(&A, b, x);
    clock_t end_time = clock();
    double cpu_time_taken = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;

    //printing the outputs
    printf("The matrix name: %s\n", filename);
    printf("The dimension of the matrix: %d by %d\n", A.num_rows, A.num_cols);
    printf("Number of non-zeros: %d\n", A.row_ptr[A.num_rows]);
    printf("CPU time taken solve Ax=b: %lf\n", cpu_time_taken);

    //computing and printing the residual and residual norm
    double *residual = (double *)malloc(A.num_rows * sizeof(double));
    compute_residual(&A, b, x, residual);
    double residual_norm = compute_norm(residual, A.num_rows);
    printf("Residual Norm: %.16e\n", residual_norm);

    // Clean up
    free(A.csr_data);
    free(A.col_ind);
    free(A.row_ptr);
    free(x);
    free(b);
    free(residual);

    return EXIT_SUCCESS;
}