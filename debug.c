#include "Matrix.h"
#include "time.h"
int main(){
    struct Matrix x, y;
    int i, j;
    time_t time1, time2;
    double cos, sin;
    int *pointer = &j;
    Matrix_scan(&x);
    //cut_last_column_and_last_string(&x);
    Matrix_print(x);
    scanf("%d", &i);
    time(&time1);
    y = Inverse_Matrix_parallel(x, i);
    time(&time2);
    printf("parallel time is %f\n", difftime(time1, time2));
    Matrix_print(y);
    scanf("%d", &i);
    time(&time1);
    y = Inverse_Matrix(x);
    time(&time2);
    printf("instead of %f usually\n", difftime(time1, time2));
    Matrix_print(y);
    //scanf("%d %d %lf %lf", &i, &j, &cos, &sin);
    //y = QR_decomposition(&x);
    //Matrix_print(y);
    //Matrix_print(x);
    //c = Matrix_multiplication(x, y);
    //Matrix_print(c);
    //c = R_Matrix_multiplication_left(x, y);
    //Matrix_print(c);
    //y = x;
    //Q = almost_triangle_form(&x);
    /*double * eigenvalues = (double* ) calloc(x.num_of_columns, sizeof(double));
    QR_algorithm(&x, eigenvalues);
    printf("%lf\n %lf\n", eigenvalues[0], eigenvalues[1]);*/
    //Matrix_print(y);
    //Matrix_print(Q);
    //Matrix_print(x);
    //x = Matrix_multiplication(Q, x);
    //Matrix_print(x);
    //free(eigenvalues);
    free(x.elements);
    //free(Q.elements);
    return 0;
}