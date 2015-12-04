#include "Matrix.h"

int main(){
    struct Matrix x, y;
    int i, j;
    double cos, sin;
    Matrix_scan(&x);
    //cut_last_column_and_last_string(&x);
    Matrix_print(x);
    //y = x;
    //Q = almost_triangle_form(&x);
    double * eigenvalues = (double* ) calloc(x.num_of_columns, sizeof(double));
    QR_algorithm(&x, eigenvalues);
    printf("%lf\n %lf\n", eigenvalues[0], eigenvalues[1]);
    //Matrix_print(y);
    //Matrix_print(Q);
    //Matrix_print(x);
    //x = Matrix_multiplication(Q, x);
    //Matrix_print(x);
    //free(eigenvalues);
    free(x.elements);
    free(y.elements);
    //free(Q.elements);
    return 0;
}