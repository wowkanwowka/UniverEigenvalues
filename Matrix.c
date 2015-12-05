#include "Matrix.h"

int main(int argc, char** argv){
    struct Matrix M;
    M.elements = NULL;
    int i;
    double m_trace, sum_of_eigenvalues;
    if((argc == 2) && !((strcmp(argv[1], "screen") != 0) && ((strcmp(argv[1], "Screen") != 0)))){
        Matrix_scan(&M);
        Matrix_print(M);
    }
    else if((argc == 3) && !((strcmp(argv[1], "file") != 0) && ((strcmp(argv[1], "File") != 0)))){
        if(Matrix_file_scan(&M, argv[2])){
            Matrix_print(M);
        }
        else{
            return 0;
        }
    }
    else if((argc == 2) && !((strcmp(argv[1], "formula") != 0) && ((strcmp(argv[1], "Formula") != 0)))){
        formula_initialization(&M);
        Matrix_print(M);
    }
    else{
        printf("Incorrect input, please enter arguments as follows:\n1) 'screen' (or 'Screen') for input from screen\n2) 'file' (or 'File') filename for input from file with name 'filename'\n3) 'formula' (or 'Formula') for formula initialization\n");
        return 0;
    }
    m_trace = 0.0;
    m_trace = trace(&M);
    
    double * eigenvalues = (double *) calloc(M.num_of_columns * M.num_of_strings, sizeof(double));
    QR_algorithm(&M, eigenvalues);
    printf("%lf\n%lf\n", eigenvalues[1], eigenvalues[0]);
    sum_of_eigenvalues = 0.0;
    for(i = 0; i < M.num_of_columns; i++){
        sum_of_eigenvalues += eigenvalues[i];
    }
    printf("nevyazka1: %lf\n", fabs(sum_of_eigenvalues - m_trace));
    m_trace = Matrix_norm_2(M);
    sum_of_eigenvalues = 0.0;
    for(i = 0; i < M.num_of_columns; i++){
        sum_of_eigenvalues += eigenvalues[i] * eigenvalues[i];
    }
    printf("nevyazka2: %lf\n", fabs(sum_of_eigenvalues - m_trace));
        /*for (i = 0; i < M.num_of_columns; i++){
            printf("%lf\n", eigenvalues[i]);
        }*/
    /*almost_triangle_form(&M);
    Matrix_print(M);
    QR_decomposition(&M);
    Matrix_print(M);*/
    /*printf("Nevyazka1 = %g\n", Matrix_norm_1(s));
    printf("Nevyazka2 = %g\n", Matrix_norm_2(s));
    printf("Nevyazka3 = %g\n", Matrix_norm_3(s));*/
    //free(s.elements);
    free(M.elements);
    return 0;
}