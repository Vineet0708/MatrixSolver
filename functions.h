#ifndef FUNCTIONS_H
#define FUNCTIONS_H

// ###########################################################
// Do not change this part
typedef struct {
    double *csr_data;   // Array of non-zero values
    int *col_ind;       // Array of column indices
    int *row_ptr;       // Array of row pointers
    int num_non_zeros;  // Number of non-zero elements
    int num_rows;       // Number of rows in matrix
    int num_cols;       // Number of columns in matrix
} CSRMatrix;


void ReadMMtoCSR(const char *filename, CSRMatrix *matrix);
void spmv_csr(const CSRMatrix *A, const double *x, double *y);

// ###########################################################

/* <Here you can add the declaration of functions you need.>
<The actual implementation must be in functions.c>
Here what "potentially" you need:

1. A function called "solver" receiving const CSRMatrix A, double *b, double *x
and solving the linear system */
void solver(const CSRMatrix *A, const double *b, double *x);
/*
2. A function called "compute_residual" to compute the residual like r=Ax-b.
This shows how much x values you found are accurate, but 
printing the whole vector r might not be a good idea. So */
void compute_residual(const CSRMatrix *A, const double *b, const double *x, double *residual);
/*
3. A function called compute_norm to compute the norm of vector residual
*/
double compute_norm(const double *vector, int size);

//A function to solve the matrix multiplication problem
void GuassSeidel(const CSRMatrix *A, const double *b, double *x, int iterations);

#endif
