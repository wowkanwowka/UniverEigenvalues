#include "Matrix.h"

void Matrix_scan(struct Matrix *x){
     int m, n, i, j;
    printf("Enter the size of matrix:");
    scanf("%d %d", &n, &m);
    while((n <= 0)|| (m <= 0)){
        printf("Are you sure?\nTry again:");
        scanf("%d %d", &n, &m);
    }
    x->elements = (double*) realloc(x->elements, sizeof(double) * n * m);
    
    printf("Now enter elements:\n");
    for (i = 0; i < n; i++){
        for(j = 0; j < m; j++){
            scanf("%lf", &(x->elements[i * m +j]));
        }
    }
    x->num_of_strings = n;
    x->num_of_columns = m;
}

bool Matrix_file_scan(struct Matrix *x, const char* filename){
     int i, j;
    int m, n;
    FILE* in = fopen(filename, "r");
    if(in == NULL){
        printf("Error: file doesn't exist\n");
        return false;
    }
    if(fscanf(in, "%d %d", &n, &m)){
        if( (n < 0) || (m < 0)){
            printf("Error: negative number of strings or columns\n");
            return false;
        }
    }
    else{
        return false;
    }
    x->elements = (double*) realloc(x->elements, sizeof(double) * n * m);
    for (i = 0; i < n; i++){
        for(j = 0; j < m; j++){
            if(fscanf(in, "%lf", &(x->elements[i * m + j]))){
                
            }
            else{
                printf("Error: not enough elements\n");
                return false;
            }
        }
    }
    x->num_of_strings = n;
    x->num_of_columns = m;
    fclose(in);
    return true;
}

void formula_initialization(struct Matrix *x){
    int i, j, n, m;
    double A = 0.0;
    double P = ((double)4)*(atan((double)1));
    printf("Enter size:");
    
    scanf("%d %d", &n, &m);
    printf("Enter A: ");
    scanf("%lf", &A);
    A *= A;
    x->elements = (double*) realloc(x->elements, sizeof(double) * n * m);
    x->num_of_columns = m;
    x->num_of_strings = n;
    for (i = 0; i < n; i++){
        for(j = 0; j < m; j++){
            if(i == j){
                // printf("%lf %lf %lf",A / 3.0, 2.0 * A  / ((2.0 * P * (double)(i + 1)) * (2.0 * P * (double)(i + 1))), (P * (double)(i + 1)) * (P * (double)(i + 1)) / A);
                x->elements[i * x->num_of_columns + j] = A / 3.0 - 2.0 * A  / ((2.0 * P * (double)(i + 1)) * (2.0 * P * (double)(i + 1))) + (P * (double)(i + 1)) * (P * (double)(i + 1)) / A;
            }
            else{
                x->elements[i * x->num_of_columns + j] = (2.0 * A  * (1.0 / ((double)(i - j) * (double)(i - j)) - 1.0 / ((double)(i + j + 2) * (double)(i + j + 2)))) / (P * P);
            }
        }
    }
}
void Matrix_print(struct Matrix x){
     int i = 0;
     int j = 0;
    for (i = 0; i < x.num_of_strings; i++){
        for(j = 0; j < x.num_of_columns; j++){
            if(fabs(x.elements[i * x.num_of_columns +j]) < EPS){
                printf("%g ", 0.0);
            }
            else{
                printf("%g ", x.elements[i * x.num_of_columns +j]);
            }
        }
        printf("\n");
    }
    printf("\n");
}

struct Matrix Matrix_addition(struct Matrix* const x, struct Matrix* const y){
    struct Matrix s;
    int i, j;
    if((x->num_of_strings != y->num_of_strings) || (x->num_of_columns != y->num_of_columns)){
        printf("Error: Addition cannot be performed!");
        return s;
    }
    s.elements = (double *) calloc(x->num_of_strings * x->num_of_columns, sizeof(double));
    s.num_of_strings = x->num_of_strings;
    s.num_of_columns = x->num_of_columns;
    for (i = 0; i < s.num_of_strings; i++){
        for(j = 0; j < s.num_of_columns; j++){
            s.elements[i * s.num_of_columns + j] = x->elements[i * x->num_of_columns + j] + y->elements[i * y->num_of_columns + j];
        }
    }
    return s;
}

struct Matrix Matrix_multiplication(struct Matrix x, struct Matrix y){
    struct Matrix s;
    int i, j, k;
    double l;
    if(x.num_of_columns != y.num_of_strings){
        printf("Error: Multiplication cannot be performed!");
        return x;
    }
    s.elements = (double*) malloc(sizeof(double) * x.num_of_strings * y.num_of_columns);
    s.num_of_strings = x.num_of_strings;
    s.num_of_columns = y.num_of_columns;
    for(i = 0; i < s.num_of_strings; i++){
        for(j = 0; j < s.num_of_columns; j++){
            l = 0;
            for(k = 0; k < y.num_of_strings; k++){
                l += x.elements[i * x.num_of_columns + k] * y.elements[k * y.num_of_columns + j];
            }
            s.elements[i * s.num_of_columns + j] = l;
        }
    }
    return s;
}

struct Matrix R_Matrix_multiplication_left(struct Matrix x, struct Matrix y){
    struct Matrix s;
    int i, j, k;
    double l;
    if(x.num_of_columns != y.num_of_strings){
        printf("Error: Multiplication cannot be performed!");
        return x;
    }
    s.elements = (double*) malloc(sizeof(double) * x.num_of_strings * y.num_of_columns);
    s.num_of_strings = x.num_of_strings;
    s.num_of_columns = y.num_of_columns;
    for(i = 0; i < s.num_of_strings; i++){
        for(j = 0; j < s.num_of_columns; j++){
            l = 0;
            for(k = i; k < y.num_of_strings; k++){
                l += x.elements[i * x.num_of_columns + k] * y.elements[k * y.num_of_columns + j];
            }
            s.elements[i * s.num_of_columns + j] = l;
        }
    }
    return s;
}


double trace(struct Matrix* x){
    int i;
    double trace = 0.0;
    for(i = 0; i < x->num_of_strings; i++){
        trace += x->elements[i * x->num_of_columns + i];
    }
    return trace;
}

void multiply_by_spin_matrix_left(struct Matrix const * const x, double cos, double sin, int i, int j){
    int l;
    double tmp_vector_1, tmp_vector_2;
    for(l = 0; l < x->num_of_strings; l++){
        tmp_vector_1 = cos * x->elements[i * x->num_of_columns + l] - sin * x->elements[j * x->num_of_columns + l];
        tmp_vector_2 = sin * x->elements[i * x->num_of_columns + l] + cos * x->elements[j * x->num_of_columns + l];
        x->elements[i * x->num_of_columns + l] = tmp_vector_1;
        x->elements[j * x->num_of_columns + l] = tmp_vector_2;
    }
}

void multiply_by_spin_matrix_right(struct Matrix const * const x, double cos, double sin, int i, int j){
    int k;
    double tmp_vector_1, tmp_vector_2;
    for(k = 0; k < x->num_of_strings; k++){
        tmp_vector_1 = cos * x->elements[k * x->num_of_columns + i] + sin * x->elements[k * x->num_of_columns + j];
        tmp_vector_2 = -sin * x->elements[k * x->num_of_columns + i] + cos * x->elements[k * x->num_of_columns + j];
        x->elements[k * x->num_of_columns + i] = tmp_vector_1;
        x->elements[k * x->num_of_columns + j] = tmp_vector_2;
    }
}

void cut_last_column_and_last_string(struct Matrix *x){
    int i, j;
    for (i = 0; i < x->num_of_strings - 1; i++){
        for(j = 0; j < x->num_of_columns - 1; j++){
            x->elements[i * (x->num_of_columns - 1) + j] = x->elements[i * x->num_of_columns + j];
        }
    }
    x->num_of_strings--;
    x->num_of_columns--;
}

void transfer_of_I_type( int i,  int j, double k, struct Matrix *x){
     int l;
    for(l=0 ; l < x->num_of_columns ; l++){
        x->elements[i * x->num_of_columns + l] -= k * x->elements[j * x->num_of_columns + l];
    }
}

void transfer_of_II_type( int i, int j, struct Matrix *x){
    double temp;
     int l;
    for(l=0; l < x->num_of_columns; l++){
        temp = x->elements[i * x->num_of_columns + l];
        x->elements[i * x->num_of_columns + l] = x->elements[j * x->num_of_columns + l];
        x->elements[j * x->num_of_columns + l]=temp;
    }
}

void transfer_of_I_type_parallel(void * argument){
    struct ARGS *  args = (struct ARGS*) argument;
    int i;
    //Matrix_print(args->matr);
    //printf("new thread %d %d %d\n", args->beginning, args->ending, args->strnum);
    for(i = args->beginning; i < args->ending; i++){
        if(fabs(args->matr.elements[args->strnum * args->matr.num_of_columns + args->strnum]) > EPS){
            transfer_of_I_type(i, args->strnum, (args->matr.elements[i * args->matr.num_of_columns + args->strnum] / args->matr.elements[args->strnum * args->matr.num_of_columns + args->strnum]), &args->matr);
        }
    }
    pthread_exit(NULL);
}

void transfer_of_III_type( int i, double k, struct Matrix *x){
     int l;
    for(l=0; l<x->num_of_columns; l++){
        x->elements[i * x->num_of_columns + l] *= k;
    }
}
void column_transfer_of_II_type( int i,  int j, struct Matrix *x){
    double temp;
     int l;
    for(l=0; l < x->num_of_strings; l++){
        temp = x->elements[l * x->num_of_columns + i];
        x->elements[l * x->num_of_columns + i] = x->elements[l * x->num_of_columns + j];
        x->elements[l * x->num_of_columns + j]=temp;
    }
}

int* Gauss_style(int *hmtt2, struct Matrix *x){
    
     int i, j, k, col_num_1, l, counter = 0;
    double max;
    int *permutations = (int *) calloc(2 * x->num_of_strings, sizeof(int));
    bool truth=true;
    *hmtt2=0;
    k = 0;
    for(i = 0; (i < x->num_of_columns) && (k < x->num_of_strings); i++){
        j = k;
        truth = true;
        max = fabs(x->elements[j * x->num_of_columns + j]);
        col_num_1 = j;
        for (l = j; l < x->num_of_strings; l++){
            if (max < fabs(x->elements[j * x->num_of_columns + l])){
                col_num_1 = l;
                max = fabs(x->elements[j * x->num_of_columns + l]);
            }
        }
        if(!(j == col_num_1)){
            column_transfer_of_II_type(j, col_num_1, x);
            permutations[2 * counter] = j + 1;
            permutations[2 * counter + 1] = col_num_1 + 1;
            counter++;
        }
        if(!(fabs(x->elements[j * x->num_of_columns + i]) < EPS)){
            truth = false;
        }
        while(truth){
            j++;
            if(!(j < x->num_of_strings)){
                truth=false;
            }
            else{
                if(!(fabs(x->elements[j * x->num_of_columns + i]) < EPS)){
                    truth=false;
                }
            }
        }
        if((j < x->num_of_strings)){
            if(j > k){
                transfer_of_II_type(k, j, x);
                (*hmtt2)++;
            }
            j++;
            for(; j < x->num_of_strings; j++){
                if(fabs(x->elements[j * x->num_of_columns + i]) > EPS){
                    transfer_of_I_type(j, k, (x->elements[j * x->num_of_columns + i] / x->elements[k * x->num_of_columns + i]), x);
                }
            }
            k++;
        }
    }
    return permutations;
}

int* Gauss_style_parallel(int *hmtt2, struct Matrix *x, int num_of_threads){
    pthread_t* threads;
    struct ARGS* threads_arguments;
    int i, j, k, col_num_1, l, counter = 0;
    double max;
    int *permutations = (int *) calloc(2 * x->num_of_strings, sizeof(int));
    bool truth=true;
    //printf("numoftreads is %d\n", num_of_threads);
    threads = (pthread_t *) malloc(sizeof(pthread_t) * num_of_threads);
    threads_arguments = (struct ARGS*) malloc(sizeof(struct ARGS) * num_of_threads);
    *hmtt2=0;
    k = 0;
    for(i = 0; (i < x->num_of_columns) && (k < x->num_of_strings); i++){
        j = k;
        truth = true;
        max = fabs(x->elements[j * x->num_of_columns + j]);
        col_num_1 = j;
        for (l = j; l < x->num_of_strings; l++){
            if (max < fabs(x->elements[j * x->num_of_columns + l])){
                col_num_1 = l;
                max = fabs(x->elements[j * x->num_of_columns + l]);
            }
        }
        if(!(j == col_num_1)){
            column_transfer_of_II_type(j, col_num_1, x);
            permutations[2 * counter] = j + 1;
            permutations[2 * counter + 1] = col_num_1 + 1;
            counter++;
        }
        if(!(fabs(x->elements[j * x->num_of_columns + i]) < EPS)){
            truth = false;
        }
        while(truth){
            j++;
            if(!(j < x->num_of_strings)){
                truth=false;
            }
            else{
                if(!(fabs(x->elements[j * x->num_of_columns + i]) < EPS)){
                    truth=false;
                }
            }
        }
        if((j < x->num_of_strings)){
            if(j > k){
                transfer_of_II_type(k, j, x);
                (*hmtt2)++;
            }
            j++;
            if((x->num_of_strings - j) / num_of_threads > 0){
                for(l = 0; l < num_of_threads; l++){
                    //printf("thread number %d\n", l);
                    threads_arguments[l].matr = *x;
                    threads_arguments[l].strnum = k;
                    threads_arguments[l].beginning = j + l * ((x->num_of_strings - j) / num_of_threads);
                    if(l < num_of_threads - 1){
                        threads_arguments[l].ending = j + (l + 1) * ((x->num_of_strings - j) / num_of_threads);
                    }
                    else{
                        threads_arguments[l].ending = x->num_of_strings;
                    }
                    pthread_create(&threads[l], NULL, transfer_of_I_type_parallel,  threads_arguments + l);
                    
                    
                }
                for(l = 0; l < num_of_threads; l++){
                    pthread_join(threads[l], NULL);
                }
                
            }
            else{
                for(; j < x->num_of_strings; j++){
                    if(fabs(x->elements[j * x->num_of_columns + i]) > EPS){
                        transfer_of_I_type(j, k, (x->elements[j * x->num_of_columns + i] / x->elements[k * x->num_of_columns + i]), x);
                    }
                }
            }
            k++;
        }
    }
    free(threads);
    free(threads_arguments);
    return permutations;
}


struct Matrix Inverse_Matrix(struct Matrix x){

    if(!(x.num_of_strings == x.num_of_columns)){
        printf("!!!SO STUPID!!!");
    }
    else{
        struct Matrix s;
        int i, j, k = 0;
        s.elements = (double*) malloc(sizeof(double) * 2 * x.num_of_strings * x.num_of_columns);
        s.num_of_strings = x.num_of_strings;
        s.num_of_columns = 2 * x.num_of_columns;
        for(i = 0; i < x.num_of_strings; i++){
            for(j=0; j < x.num_of_columns; j++){
                s.elements[i * 2 * x.num_of_columns + j] = x.elements[i * x.num_of_columns + j];
            }
            s.elements[i * s.num_of_columns + x.num_of_columns + i] = 1.0;
        }
        int* transfers = Gauss_style(&k, &s);
        for(i = 0; i < x.num_of_strings; i++){
            if (fabs(s.elements[i * s.num_of_columns + i]) < EPS){
                printf("%d %lf", i, s.elements[i * s.num_of_columns + i]);
                free(s.elements);
                printf("Error matrix is degenerate argument will be returned\n");
                return x;
            }
        }
        transfer_of_III_type(x.num_of_strings - 1, 1.0 / (s.elements[(x.num_of_strings-1) * s.num_of_columns + x.num_of_strings - 1]), &s);
        for(i = x.num_of_strings - 1; i > -1; i--){
            for(j = i - 1; j > -1; j--){
                if(fabs(s.elements[j * s.num_of_columns + i]) > 0){
                    transfer_of_I_type(j, i, s.elements[j * s.num_of_columns + i] / s.elements[i * s.num_of_columns + i], &s);
                }
            }
        }
        for(i = 0; i < x.num_of_strings - 1; i++){
            transfer_of_III_type(i, 1.0 / (s.elements[i * s.num_of_columns + i]), &s);
        }
        for(i = 0; (i < x.num_of_strings) && (transfers[2 * i] > 0); i++){
        }
        for(i -= 1; (i >= 0) && (transfers[2 * i] > 0); i--){
            transfer_of_II_type(transfers[2 * i] - 1, transfers[2 * i + 1] - 1, &s);
        }
        struct Matrix w;
        w.elements = (double*) malloc(sizeof(double) * x.num_of_strings * x.num_of_columns);
        w.num_of_columns = x.num_of_columns;
        w.num_of_strings = x.num_of_strings;
        for(i = 0; i < x.num_of_strings; i++){
            for(j = 0; j < x.num_of_strings; j++){
                w.elements[i * w.num_of_columns + j] = s.elements[i * s.num_of_columns + x.num_of_strings + j];
            }
        }
        free(s.elements);
        return w;
    }
    return x;
}

struct Matrix Inverse_Matrix_parallel(struct Matrix x, int num_of_threads){
    pthread_t* threads;
    int l;
    struct Matrix s;
    int i, j, k = 0;
    struct ARGS* threads_arguments;
    //printf("numoftreads is %d\n", num_of_threads);
    threads = (pthread_t *) malloc(sizeof(pthread_t) * num_of_threads);
    threads_arguments = (struct ARGS*) malloc(sizeof(struct ARGS) * num_of_threads);

    if(!(x.num_of_strings == x.num_of_columns)){
        printf("!!!SO STUPID!!!");
    }
    else{
        
        s.elements = (double*) malloc(sizeof(double) * 2 * x.num_of_strings * x.num_of_columns);
        s.num_of_strings = x.num_of_strings;
        s.num_of_columns = 2 * x.num_of_columns;
        for(i = 0; i < x.num_of_strings; i++){
            for(j=0; j < x.num_of_columns; j++){
                s.elements[i * 2 * x.num_of_columns + j] = x.elements[i * x.num_of_columns + j];
            }
            s.elements[i * s.num_of_columns + x.num_of_columns + i] = 1.0;
        }
        int* transfers = Gauss_style_parallel(&k, &s, num_of_threads);
        for(i = 0; i < x.num_of_strings; i++){
            if (fabs(s.elements[i * s.num_of_columns + i]) < EPS){
                printf("%d %lf", i, s.elements[i * s.num_of_columns + i]);
                free(s.elements);
                printf("Error matrix is degenerate argument will be returned\n");
                return x;
            }
        }
        transfer_of_III_type(x.num_of_strings - 1, 1.0 / (s.elements[(x.num_of_strings-1) * s.num_of_columns + x.num_of_strings - 1]), &s);
        for(i = x.num_of_strings - 1; i > -1; i--){
            if((i + 1) / num_of_threads > 0){
                for(l = 0; l < num_of_threads; l++){
                    //printf("thread number %d\n", l);
                    threads_arguments[l].matr = s;
                    threads_arguments[l].strnum = i;
                    threads_arguments[l].beginning =  l * ((i + 1) / num_of_threads);
                    if(l < num_of_threads - 1){
                        threads_arguments[l].ending = (l + 1) * ((i + 1) / num_of_threads);
                    }
                    else{
                        threads_arguments[l].ending = i;
                    }
                    pthread_create(&threads[l], NULL, transfer_of_I_type_parallel,  threads_arguments + l);
                    
                    
                }
                for(l = 0; l < num_of_threads; l++){
                    pthread_join(threads[l], NULL);
                }
                
            }
            else{
                for(j = i - 1; j > -1; j--){
                    if(fabs(s.elements[j * s.num_of_columns + i]) > 0){
                        transfer_of_I_type(j, i, s.elements[j * s.num_of_columns + i] / s.elements[i * s.num_of_columns + i], &s);
                    }
                }
            }
        }
        free(threads_arguments);
        free(threads);
        for(i = 0; i < x.num_of_strings - 1; i++){
            transfer_of_III_type(i, 1.0 / (s.elements[i * s.num_of_columns + i]), &s);
        }
        for(i = 0; (i < x.num_of_strings) && (transfers[2 * i] > 0); i++){
        }
        for(i -= 1; (i >= 0) && (transfers[2 * i] > 0); i--){
            transfer_of_II_type(transfers[2 * i] - 1, transfers[2 * i + 1] - 1, &s);
        }
        struct Matrix w;
        w.elements = (double*) malloc(sizeof(double) * x.num_of_strings * x.num_of_columns);
        w.num_of_columns = x.num_of_columns;
        w.num_of_strings = x.num_of_strings;
        for(i = 0; i < x.num_of_strings; i++){
            for(j = 0; j < x.num_of_strings; j++){
                w.elements[i * w.num_of_columns + j] = s.elements[i * s.num_of_columns + x.num_of_strings + j];
            }
        }
        free(s.elements);
        return w;
    }
    return x;
}


double Matrix_norm_1(struct Matrix x){
     int i,j;
    double norm = 0.0;
    for(i = 0; i < x.num_of_strings; i++){
        for(j = 0; j < x.num_of_columns; j++){
            norm += fabs(x.elements[i * x.num_of_columns + j]);
        }
    }
    return norm;
}

double Matrix_norm_2(struct Matrix x){
     int i,j;
    double norm = 0.0;
    for(i = 0; i < x.num_of_strings; i++){
        for(j = 0; j < x.num_of_columns; j++){
            norm += x.elements[i * x.num_of_columns + j] * x.elements[i * x.num_of_columns + j];
        }
    }
    return norm;
}

double Matrix_norm_3(struct Matrix x){
     int i,j;
    double norm = 0.0;
    for(i = 0; i < x.num_of_strings; i++){
        for(j = 0; j < x.num_of_columns; j++){
            if(fabs(x.elements[i * x.num_of_columns + j]) > norm){
                norm = fabs(x.elements[i * x.num_of_columns + j]);
            }
        }
    }
    return norm;
}

struct Matrix QR_decomposition(struct Matrix *x){
    //printf("decomposing\n");
    struct Matrix Q, TMP;
    int i, j, k;
    double tmp, cos, sin;
    //Matrix_print(x);
    Q.num_of_columns = x->num_of_columns;
    Q.num_of_strings = x->num_of_strings;
    Q.elements = (double *) calloc(x->num_of_strings * x->num_of_columns, sizeof(double));
    for(k = 0; k < Q.num_of_strings; k++){
        Q.elements[k * Q.num_of_columns + k] = 1.0;
    }
    for(i = 0; i < x->num_of_strings - 1; i++){
        for(j = i + 1; (j < x->num_of_strings) && (fabs(x->elements[j * x->num_of_columns + i]) > EPS); j++){
            if((fabs(x->elements[i * x->num_of_columns + i]) < EPS) && (fabs(x->elements[j * x->num_of_columns + i]) < EPS)){
                cos = 1.0;
                sin = 0.0;
            }
            else{
                tmp = sqrt(x->elements[i * x->num_of_columns + i] * x->elements[i * x->num_of_columns + i] + x->elements[j * x->num_of_columns + i] * x->elements[j * x->num_of_columns + i]);
                cos = x->elements[i * x->num_of_columns + i] / tmp;
                
                sin = -x->elements[j * x->num_of_columns + i] / tmp;
                multiply_by_spin_matrix_left(&Q, cos, sin, i, j);
                multiply_by_spin_matrix_left(x, cos, sin, i, j);
                //printf("%lf %lf\n", cos, sin);
            }
            
        }
        
    }
    for(i = 0; i < Q.num_of_strings; i++){
        for(j = i + 1; j < Q.num_of_columns; j++){
            tmp = Q.elements[i * Q.num_of_columns + j];
            Q.elements[i * Q.num_of_columns + j] = Q.elements[j * Q.num_of_columns + i];
            Q.elements[j * Q.num_of_columns + i] = tmp;
        }
    }
    return Q;
}

struct Matrix QR_decomposition_reflections(struct Matrix *x){
    struct Matrix s;
    struct Matrix Q;
    struct Matrix tmp_vector;
    int i, j, k;
    double tmp = 0.0;
    s.num_of_columns = x->num_of_columns;
    s.num_of_strings = x->num_of_strings;
    Q.num_of_columns = x->num_of_columns;
    Q.num_of_strings = x->num_of_strings;
    Q.elements = (double *) calloc(x->num_of_strings * x->num_of_columns, sizeof(double));
    for(k = 0; k < Q.num_of_strings; k++){
        Q.elements[k * Q.num_of_columns + k] = 1.0;
    }
    tmp_vector.elements = (double*) calloc(x->num_of_strings, sizeof(double));
    tmp_vector.num_of_strings = x->num_of_strings;
    tmp_vector.num_of_columns = 1;
    s.elements = (double*) calloc(x->num_of_strings, sizeof(double));
    for(i = 0; i < x->num_of_strings - 1; i++){
        tmp = 0.0;
        for(j = i; j < x->num_of_strings; j++){
            tmp_vector.elements[j] = x->elements[j * x->num_of_columns + i];
            if((fabs(x->elements[j * x->num_of_columns + i]) > tmp) && (j > i)){
                tmp = fabs(x->elements[j * x->num_of_columns + i]);
            }
        }
        if(tmp > EPS){
            //printf("%d ", i);
            //Matrix_print(tmp_vector);
            s = QR_reflection(&tmp_vector, i);
            //Matrix_print(s);
            //printf("%d %d\n", i, j);
            //Matrix_print(s);
            //*x = Matrix_multiplication(s, *x);
            Q = Matrix_multiplication(s, Q);
            *x = Matrix_multiplication(s, *x);
            //Matrix_print(*x);
            //Matrix_print(s);
        }
    }
    free(s.elements);
    //s = Matrix_multiplication(s, x);
    //Matrix_print(s);
    //Matrix_print(*x);
    //Matrix_print(Q);
    for(i = 0; i < Q.num_of_strings; i++){
        for(j = i + 1; j < Q.num_of_columns; j++){
            tmp = Q.elements[i * Q.num_of_columns + j];
            Q.elements[i * Q.num_of_columns + j] = Q.elements[j * Q.num_of_columns + i];
            Q.elements[j * Q.num_of_columns + i] = tmp;
        }
    }
    return Q;
}

struct Matrix Multiply_by_number(struct Matrix* const x, double k){
    struct Matrix s;
    int i, j;
    s.elements = (double*) calloc(x->num_of_strings * x->num_of_columns, sizeof(double));
    s.num_of_columns = x->num_of_columns;
    s.num_of_strings = x->num_of_strings;
    for(i = 0; i < s.num_of_strings; i++){
        for(j = 0; j < s.num_of_columns; j++){
            s.elements[i * s.num_of_columns + j]=k * x->elements[i * x->num_of_columns + j];
        }
    }
    return s;
}


struct Matrix almost_triangle_form(struct Matrix *x){
    struct Matrix Q;
    int i, j, k, counter;
    counter = 0;
    double tmp, cos, sin;
    //Matrix_print(x);
    /*Q.num_of_columns = x->num_of_columns;
    Q.num_of_strings = x->num_of_strings;
    Q.elements = (double *) calloc(x->num_of_strings * x->num_of_columns, sizeof(double));
    for(k = 0; k < Q.num_of_strings; k++){
        Q.elements[k * Q.num_of_columns + k] = 1.0;
    }*/
    for(i = 1; i < x->num_of_strings - 1; i++){
        for(j = i + 1; j < x->num_of_strings; j++){
                tmp = sqrt(x->elements[i * x->num_of_columns + i - 1] * x->elements[i * x->num_of_columns + i - 1] + x->elements[j * x->num_of_columns + i - 1] * x->elements[j * x->num_of_columns + i - 1]);
                if(tmp > EPS){
                cos = x->elements[i * x->num_of_columns + i - 1] / tmp;
                sin = -x->elements[j * x->num_of_columns + i - 1] / tmp;
                    if((fabs(sin) > EPS)){
                        //printf("step: %d %d\n", i, j);
                        //printf("spin_mult %d\n", counter++);
                        multiply_by_spin_matrix_left(x, cos, sin, i, j);
                        //Q = multiply_by_spin_matrix_left(&Q, cos, sin, i, j);
                        //printf("spin_mult %d\n", counter++);
                        multiply_by_spin_matrix_right(x, cos, -sin, i, j);
                    }
            }
        }
        
    }
    /*for(i = 0; i < Q.num_of_strings; i++){
        for(j = i + 1; j < Q.num_of_columns; j++){
            tmp = Q.elements[i * Q.num_of_columns + j];
            Q.elements[i * Q.num_of_columns + j] = Q.elements[j * Q.num_of_columns + i];
            Q.elements[j * Q.num_of_columns + i] = tmp;
        }
    }*/
    return Q;
}

void QR_algorithm(struct Matrix * const x, double* eigenvalues){
    struct Matrix A, Q, TMP;
    int counter = 0;
    int k = 0;
    int i, j;
    double sdvig = 0.0;
    double a, b, c;
    bool flag = true;
    if(x->num_of_columns == 2){
        c = x->elements[0] * x->elements[3] - x->elements[1] * x->elements[2];
        b = x->elements[0] + x->elements[3];
        //printf("%lf %lf\n", b, c);
        if(b * b - 4 * c > 0){
            c = sqrt(b * b - 4 * c);
            eigenvalues[0] = (b + c) / 2.0;
            eigenvalues[1] = (b - c) / 2.0;
        }
        else if(fabs(b * b - 4 * c) < EPS){
            eigenvalues[0] = b / 2.0;
            eigenvalues[1] = b / 2.0;
        }
        else{
            printf("complex eigenvalues\n");
            return ;
        }
    }
    else{
        A.elements = (double*) calloc(x->num_of_strings * x->num_of_columns, sizeof(double));
        A.num_of_columns = x->num_of_columns;
        A.num_of_strings = x->num_of_strings;
        for(i = 0; i < A.num_of_strings; i++){
            for(j = 0; j  < A.num_of_columns; j++){
                A.elements[i * A.num_of_columns + j] = x->elements[i * x->num_of_columns + j];
            }
        }
        almost_triangle_form(&A);
        //Matrix_print(A);
        for(i = 0; i < x->num_of_columns - 2; i++){
            flag = true;
            while(flag){
                
                if(fabs(A.elements[A.num_of_columns * A.num_of_strings - 2]) < EPS){
                    //Matrix_print(A);
                    
                    
                    //Matrix_print(TMP);
                    printf("%lf\n", A.elements[A.num_of_strings * A.num_of_columns - 1]);
                    eigenvalues[A.num_of_strings - 1] = A.elements[A.num_of_strings * A.num_of_columns - 1];
                    cut_last_column_and_last_string(&A);
                    //printf("%d %d %d\n", A.num_of_columns, A.num_of_strings, i);
                    flag = false;
                    //printf("Do you wish to proceed?\n");
                    //scanf("%d", &j);
                }
                if(flag){
                    //printf("new iteration step:%d\n", i);
                    if(counter > 2000){
                        printf("Too many iterations required, here's what we've got");
                        printf("%lf\n", A.elements[A.num_of_strings * A.num_of_columns - 1]);
                        eigenvalues[A.num_of_columns - 1] = A.elements[A.num_of_columns * A.num_of_strings - 1];
                        cut_last_column_and_last_string(&A);
                        flag = false;
                    }
                    else{
                        if(fabs(A.elements[A.num_of_strings * A.num_of_columns - 1]) < EPS){
                            sdvig = A.elements[A.num_of_strings * A.num_of_columns - 2];
                        }
                        else{
                            sdvig = A.elements[A.num_of_strings * A.num_of_columns - 1];
                        }
                        //Matrix_print(A);
                        for(k = 0; k < A.num_of_strings; k++){
                            A.elements[k * A.num_of_columns + k] -= sdvig;
                        }
                        //Matrix_print(A);
                        Q = QR_decomposition(&A);
                        //Matrix_print(A);
                        //Matrix_print(Q);
                        /*for(i = 0; i < Q.num_of_strings; i++){
                            for(j = 0; j < Q.num_of_columns; j++){
                                Q1.elements[i * Q.num_of_columns + j] = Q.elements[j * Q.num_of_columns + i];
                            }
                        }*/
                        TMP = R_Matrix_multiplication_left(A, Q);
                        free(A.elements);
                        A = TMP;
                        for(k = 0; k < A.num_of_strings; k++){
                            A.elements[k * A.num_of_columns + k] += sdvig;
                        }
                        //Matrix_print(A);
                        //Matrix_print(A);
                        free(Q.elements);
                        counter++;
                    }
                }
            }
        }
        c = A.elements[0] * A.elements[3] - A.elements[1] * A.elements[2];
        b = A.elements[0] + A.elements[3];
        //printf("%lf %lf\n", b, c);
        if(b * b - 4 * c > 0){
            c = sqrt(b * b - 4 * c);
            eigenvalues[0] = (b + c) / 2.0;
            eigenvalues[1] = (b - c) / 2.0;
        }
        else if(fabs(b * b - 4 * c) < EPS){
            eigenvalues[0] = b / 2.0;
            eigenvalues[1] = b / 2.0;
        }
        else{
            printf("complex eigenvalues\n");
            return ;
        }

    }
    free(A.elements);
}

struct Matrix reflection_matrix(struct Matrix const * const x){
    int i;
    struct Matrix E;
    struct Matrix result, x_transposed;
    result.num_of_columns = 0;
    result.num_of_strings = 0;
    if(x->num_of_columns != 1){
        result.elements = NULL;
        printf("Error: argument is not a vector\n");
        return result;
    }
    x_transposed.num_of_strings = 1;
    x_transposed.num_of_columns = x->num_of_strings;
    x_transposed.elements = (double*) calloc(x->num_of_strings, sizeof(double));
    for(i = 0; i < x->num_of_strings; i++){
        x_transposed.elements[i] = x->elements[i];
    }
    E.elements = (double *) calloc(x->num_of_strings * x->num_of_strings, sizeof(double));
    E.num_of_columns = x->num_of_strings;
    E.num_of_strings = x->num_of_strings;
    for(i = 0; i < E.num_of_strings; i++){
        E.elements[i * E.num_of_columns + i] = 1.0;
    }
    result = Matrix_multiplication(*x, x_transposed);
    result = Multiply_by_number(&result, -2.0);
    result = Matrix_addition(&E, &result);
    return result;
}

struct Matrix QR_reflection(struct Matrix const * const x, int i){
    struct Matrix tmp_vector, result;
    int j;
    if(x->num_of_columns != 1){
        result.elements = NULL;
        printf("Error: argument is not a vector\n");
        return result;
    }
    tmp_vector.num_of_columns = x->num_of_columns;
    tmp_vector.num_of_strings = x->num_of_strings;
    tmp_vector.elements = (double*) calloc(x->num_of_strings, sizeof(double));
    for(j = i; j < x->num_of_strings; j++){
        tmp_vector.elements[j] = x->elements[j];
    }
    tmp_vector.elements[i] -= Matrix_norm_2(tmp_vector);
    tmp_vector = Multiply_by_number(&tmp_vector, 1.0 / Matrix_norm_2(tmp_vector));
    //Matrix_print(tmp_vector);
    result = reflection_matrix(&tmp_vector);
    free(tmp_vector.elements);
    return result;
}

double* QR_algorithm_reflections(struct Matrix * const x){
    struct Matrix A, Q, TMP, E;
    int counter = 0;
    int k = 0;
    int i, j;
    double epsilonizator = 0.0;
    double sdvig = 0.0;
    bool flag = true;
    double* eigenvalues = (double *) calloc(x->num_of_strings, sizeof(double));
    E.elements = (double *) calloc(x->num_of_strings * x->num_of_columns, sizeof(double));
    E.num_of_columns = x->num_of_columns;
    E.num_of_strings = x->num_of_strings;
    Q.num_of_columns = x->num_of_columns;
    Q.num_of_strings = x->num_of_strings;
    A = *x;
    //Matrix_print(A);
    Q.elements = (double *) calloc(x->num_of_strings * x->num_of_columns, sizeof(double));
    for(k = 0; k < Q.num_of_strings; k++){
        Q.elements[k * Q.num_of_columns + k] = 1.0;
        E.elements[k * E.num_of_columns + k] = 1.0;
    }
    //Matrix_print(E);
    almost_triangle_form(&A);
    /*for(i = 0; i < Q.num_of_strings; i++){
     for(j = 0; j < Q.num_of_columns; j++){
     Q1.elements[i * Q.num_of_columns + j] = Q.elements[j * Q.num_of_columns + i];
     }
     }
     Matrix_print(Matrix_multiplication(Matrix_multiplication(Q, A), Q1));*/
    //Matrix_print(A);
    while(flag){
        //printf("counter is %d\n", counter);
        if(counter > 20000){
            printf("Too many iterations required, here's what we've got:");
            for(k = 0; k < x->num_of_columns; k++){
                eigenvalues[k] = A.elements[k * A.num_of_columns + k];
            }
            return eigenvalues;
        }
        if(fabs(A.elements[A.num_of_strings * A.num_of_columns - 1]) < EPS){
            sdvig = A.elements[A.num_of_strings * A.num_of_columns - 2];
        }
        else{
            sdvig = A.elements[A.num_of_strings * A.num_of_columns - 1];
        }
        //Matrix_print(A);
        TMP = Multiply_by_number(&E, -sdvig);
        //Matrix_print(TMP);
        A = Matrix_addition(&A, &TMP);
        //Matrix_print(A);
        Q = QR_decomposition_reflections(&A);
        //Matrix_print(A);
        /*for(i = 0; i < Q.num_of_strings; i++){
         for(j = 0; j < Q.num_of_columns; j++){
         Q1.elements[i * Q.num_of_columns + j] = Q.elements[j * Q.num_of_columns + i];
         }
         }*/
        A = Matrix_multiplication(A, Q);
        //Matrix_print(A);
        TMP = Multiply_by_number(&E, sdvig);
        //Matrix_print(TMP);
        A = Matrix_addition(&A, &TMP);
        //Matrix_print(A);
        //Matrix_print(A);
        for(k = 0; k < A.num_of_columns - 1; k++){
            if(fabs(A.elements[(k + 1) * A.num_of_columns + k]) > epsilonizator){
                epsilonizator = fabs(A.elements[(k + 1) * A.num_of_columns + k]);
            }
        }
        if(epsilonizator < EPS){
            flag = false;
        }
        epsilonizator = 0.0;
        counter++;
    }
    for(k = 0; k < x->num_of_columns; k++){
        eigenvalues[k] = A.elements[k * A.num_of_columns + k];
    }
    return eigenvalues;
}

