#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "../../OneDrive/Desktop/2MP3 A3/functions.h"

void ReadMMtoCSR(const char *filename, CSRMatrix *matrix) {

    // opening the file
    FILE *file = fopen(filename, "r");
    if (!file)
    {
        fprintf(stderr, "Error opening file %s\n", filename);
        exit(EXIT_FAILURE);
    }

    // reading the header line and ignoring the lines of comments
    char line[1024];
    fgets(line, sizeof(line), file);
    while (line[0] == '%') {
        fgets(line, sizeof(line), file);
    }

    // getting the dimensions from the header line
    sscanf(line, "%d %d %d", &matrix->num_rows, &matrix->num_cols, &matrix->num_non_zeros);

    //total number of non-zeros in the matrix since only symmetric part is included
    int nnz_total = 2 * (matrix->num_non_zeros) - (matrix->num_cols);

    // Allocating the memory for the CSR matrix
    matrix->csr_data = (double *) malloc(nnz_total * sizeof(double));
    matrix->col_ind = (int *) malloc(nnz_total * sizeof(int));
    matrix->row_ptr = (int *) malloc((nnz_total) * sizeof(int));

    int row, col;
    double value;

    // reading the data from the file
    for (int i = 0; i < nnz_total; i++)
    {
        fscanf(file, "%d %d %lf", &row, &col, &value);

        //assigning the vlaeus to the correct arrays
        matrix->csr_data[i] = value;
        matrix->col_ind[i] = col - 1;
        matrix->row_ptr[i] = row - 1;

        //adding the symmetric parts of the matrix into the correct columns and rows
        if (matrix->col_ind[i] != matrix->row_ptr[i])
        {
            matrix->col_ind[i+1] = matrix->row_ptr[i];
            matrix->row_ptr[i+1] = matrix->col_ind[i];
            matrix->csr_data[i+1] = matrix->csr_data[i];
            i++;
        }
    }

    int temp, min_index;
    double temp_csr;

    //bubble sort to sort in terms of rows first, and then columns
    for (int i = 0; i < nnz_total - 1; i++) {
        for (int j = 0; j < nnz_total - i - 1; j++) {
            if (matrix->row_ptr[j] > matrix->row_ptr[j + 1]) {
                // Swap rows, columns, and indices
                temp = matrix->row_ptr[j];
                matrix->row_ptr[j] = matrix->row_ptr[j + 1];
                matrix->row_ptr[j + 1] = temp;

                temp = matrix->col_ind[j];
                matrix->col_ind[j] = matrix->col_ind[j + 1];
                matrix->col_ind[j + 1] = temp;

                temp_csr = matrix->csr_data[j];
                matrix->csr_data[j] = matrix->csr_data[j + 1];
                matrix->csr_data[j + 1] = temp_csr;
            }
        }
    }

    // making the row ptrs list like it is supposed to be, not the row of the value
    int index = 0;
    matrix->row_ptr[index] = 0;
    for (int i = 1; i <= nnz_total; i++) {
        if (matrix->row_ptr[i] != matrix->row_ptr[i - 1]) {
            index++;
            matrix->row_ptr[index] = i;
        }
    }

    // reallocating the memory, making it smaller than before since the list is smaller
    index++;
    matrix->row_ptr = (int *)realloc(matrix->row_ptr, (index + 1) * sizeof(int));
    matrix->row_ptr[index] = nnz_total;

    // closing the file
    fclose(file);
}

void spmv_csr(const CSRMatrix *A, const double *x, double *y)
{
    // Iterate through each row of the matrix
    for (int i = 0; i < A->num_rows; ++i)
    {
        // Initialize the result for the current row
        double sum = 0.0;

        // Iterate through the non-zero elements of the current row
        for (int j = A->row_ptr[i]; j < A->row_ptr[i + 1]; ++j)
        {
            // Perform the multiplication and accumulate the result
            int col_index = A->col_ind[j];

            // For the lower triangular part, only consider elements where col_index <= i
            if (col_index <= i)
            {
                sum += A->csr_data[j] * x[col_index];
            }
            else
            {
                // For the upper triangular part, use the symmetry property
                sum += A->csr_data[j] * x[i];
            }
        }

        // Store the result in the output vector
        y[i] = sum;
    }
}

void GuassSeidel(const CSRMatrix *A, const double *b, double *x, int iterations) {
    // Perform Gauss-Seidel iteration for a specified number of iterations
    for (int k = 0; k < iterations; k++) {
        // Iterate through each row of the matrix
        for (int i = 0; i < A->num_rows; i++) {
            double sum = 0.0, diagonal = 0.0;

            // Calculate the sum of the non-diagonal elements in the current row
            for (int j = A->row_ptr[i]; j < A->row_ptr[i + 1]; j++) {
                if (i == A->col_ind[j]) {
                    diagonal = A->csr_data[j];
                } else {
                    sum += A->csr_data[j] * x[A->col_ind[j]];
                }
            }

            // Check if the diagonal element is close to zero
            if (diagonal == 0.0) {
                fprintf(stderr, "Error: Diagonal element is zero.\n");
                return;
            }

            // Update the solution for the current row using Gauss-Seidel iteration
            x[i] = (b[i] - sum) / diagonal;
        }
    }
}

void solver(const CSRMatrix *A, const double *b, double *x)
{
    // Set the maximum number of iterations
    int max_iterations = 1000000;

    // Call the Guass Seidel solver
    GuassSeidel(A, b, x, max_iterations);
}

void compute_residual(const CSRMatrix *A, const double *b, const double *x, double *residual) {
    // Compute the residual: residual = b - Ax
    for (int i = 0; i < A->num_rows; i++) {
        residual[i] = b[i];

        // Iterate through the non-zero elements of the current row
        for (int j = A->row_ptr[i]; j < A->row_ptr[i + 1]; j++) {
            int col = A->col_ind[j];
            residual[i] -= A->csr_data[j] * x[col];
        }
    }
}

double compute_norm(const double *vector, int size)
{
    // Compute the Euclidean norm
    double sum_of_squares = 0.0;
    for (int i = 0; i < size; ++i)
    {
        sum_of_squares += vector[i] * vector[i];
    }
    return sqrt(sum_of_squares);
}
