#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <pthread.h>
#include <time.h>

#define EPS 6e-12

struct Matrix{
    double *elements;
     int num_of_strings;
     int num_of_columns;
    
};

struct ARGS{
    struct Matrix  matr;
    int strnum, beginning, ending;
};

void Matrix_scan(struct Matrix *x);//scans Matrix from standart input stream

bool Matrix_file_scan(struct Matrix *x, const char* filename);//scans Matrix from file with name filename

void formula_initialization(struct Matrix *x);//

void Matrix_print(struct Matrix x);//prints Matrix

//struct Matrix Matrix_addition(struct Matrix x, struct Matrix y);//adds two Matrix

struct Matrix Matrix_multiplication(struct Matrix x, struct Matrix y);//multiplies two Matrix

struct Matrix R_Matrix_multiplication_left(struct Matrix x, struct Matrix y);

double trace(struct Matrix * x);

void multiply_by_spin_matrix_left(struct Matrix const * const x, double cos, double sin, int i, int j);

void multiply_by_spin_matrix_right(struct Matrix const * const x, double cos, double sin, int i, int j);

void cut_last_column_and_last_string(struct Matrix *x);

void transfer_of_I_type( int i,  int j, double k, struct Matrix *x);//subtracts j-string multiplied by k from i-string

void transfer_of_II_type( int i,  int j, struct Matrix *x);//exchanges i and j strings

void transfer_of_I_type_parallel(void * ARGS);

void column_transfer_of_II_type( int i,  int j, struct Matrix *x);//exchanges i and j columns

void transfer_of_III_type( int i, double k, struct Matrix *x);//multiplies i string by k

int* Gauss_style(int *k, struct Matrix *x);//makes Matrix Gauss styled

int* Gauss_style_parallel(int *k, struct Matrix*x, int num_of_threads);

double Determinant(struct Matrix x);//returns determinant

struct Matrix Inverse_Matrix(struct Matrix x);//returns Inverse Matrix

struct Matrix Inverse_Matrix_parallel(struct Matrix x, int num_of_threads);

struct Matrix Transpose(struct Matrix *x);//returns transpose of Matrix

//struct Matrix Multiply_by_number(double k, struct Matrix *x);//multiplies Matrix by number

double Matrix_norm_1(struct Matrix x);//returns Matrix norm (sum of absolute values)

double Matrix_norm_2(struct Matrix x);//returns Matrix norm (sqrt of sum of squares)

double Matrix_norm_3(struct Matrix x);//returns Matrix norm (max absolute value)

struct Matrix QR_decomposition(struct Matrix *x);//Returns Q matrix from QR-decomposition(spins) and change the input Matrix x

struct Matrix almost_triangle_form(struct Matrix *x);//Retruns Q matrix of basis change and makes x almost triangle

void QR_algorithm(struct Matrix * const x, double* eigenvalues);//returns pointer to the list of eigenvalues of x (stable work guaranteed only if eigenvalues are real)

double* QR_algorithm_reflections(struct Matrix * const x);//returns pointer to the list of eigenvalues of x (stable work guaranteed only if eigenvalues are real)

struct Matrix Multiply_by_number(struct Matrix* const x, double k);//Multiplies Matrix by number

struct Matrix Matrix_addition(struct Matrix* const x, struct Matrix* const y);//Adds two Matrix

struct Matrix QR_reflection(struct Matrix const * const x, int i);// Returns Q matrix from QR-decomposition(reflections) and change the input Matrix x

struct Matrix reflection_matrix(struct Matrix const * const x);//returns reflection matrix with respect to x (supposed to be of n * 1 size)
