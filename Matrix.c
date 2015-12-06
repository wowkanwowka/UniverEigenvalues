#include "Matrix.h"

int main(int argc, char** argv){
    struct Matrix M, M1, s;
    M.elements = NULL;
    int i;
    time_t time1, time2;
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
    printf("Enter number of threads: ");
    scanf("%d", &i);
    time(&time1);
    s = Inverse_Matrix_parallel(M, i);
    time(&time2);
    
    Matrix_print(s);
    s = Matrix_multiplication(s, M);
    for (i = 0; i < M.num_of_columns; i++){
        s.elements[i * s.num_of_columns + i] -= 1.0;
        }
    /*almost_triangle_form(&M);
    Matrix_print(M);
    QR_decomposition(&M);
    Matrix_print(M);*/
    printf("parallel time is %f\n", difftime(time2, time1));
    printf("Nevyazka1 = %g\n", Matrix_norm_1(s));
    printf("Nevyazka2 = %g\n", Matrix_norm_2(s));
    printf("Nevyazka3 = %g\n", Matrix_norm_3(s));
    free(s.elements);
    free(M.elements);
    return 0;
}