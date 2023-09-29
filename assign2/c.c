#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// Allocation
//////////////////////////////////////////////////////////////////
float* allocate1DArray(int n){
    float* array = (float*)malloc(n * sizeof(float));
    return array;
}

float** allocate2DArray(int rows, int cols) {
    /*
    float** a = allocate2DArray(r,c);
    */
    float** matrix = (float**)malloc(rows * sizeof(float*));
    if (matrix == NULL) {
        perror("Memory allocation failed");
        exit(1); 
    }
    for (int i = 0; i < rows; i++) {
        matrix[i] = (float*)malloc(cols * sizeof(float));
        if (matrix[i] == NULL) {
            perror("Memory allocation failed");
            exit(1); 
        }
        for (int j = 0; j < cols; j++){
            matrix[i][j] = 0.0;
        }   
    }
    return matrix;
}

// Nothing
//////////////////////////////////////////////////////////////////
void print1DMatrix(float *a, int n){
    for (int i = 0; i < n; i++) { 
        printf("%f\n", a[i]);

    }
}

void print2DMatrix(float **a, int n){
    for (int i = 0; i < n; i++) { 
        for (int j = 0; j < n; j++) {
            printf("%f\t", a[i][j]);
        }
        printf("\n"); 
    }
}

// Solve the system of equations
//////////////////////////////////////////////////////////////////
void luDecomp(float** a, float** l, float** u, int n) {
    for (int i = 0; i < n; i++){
        l[i][i] = 1.0;
    }
    for (int k = 0; k < n; k++) {

        for (int j = k; j < n; j++) {
            float sum = 0.0;
            for (int m = 0; m < k; m++) {
                sum += l[k][m] * u[m][j];
            }
            u[k][j] = a[k][j] - sum;
        }
        for (int i = k + 1; i < n; i++) {
            float sum = 0.0;
            for (int m = 0; m < k; m++) {
                sum += l[i][m] * u[m][k];
            }
            l[i][k] = (1.0 / u[k][k]) * (a[i][k] - sum);
        }
    }
}

void forwardSubstitution(float** l, float* b, float* y, int n) {
    
    y[0] = b[0] / l[0][0];
    for (int i = 1; i < n; i++) {
        float sum = 0.0;
        for (int j = 0; j < i; j++) {
            sum += l[i][j] * y[j];
        }
        y[i] = (1.0 / l[i][i]) * (b[i] - sum);
    }
}

void backwardSubstitution(float** u, float* y, float* x, int n) {
    x[n-1] = y[n-1]/u[n-1][n-1];
    for (int i = n - 2; i >= 0; i--) {
        float sum = 0.0;
        for (int j = i + 1; j < n; j++) {
            sum += u[i][j] * x[j];
        }
        x[i] = (1.0 / u[i][i]) * (y[i] - sum);
    }
}

void linalgsolve(float** a, float* b, int n){
    float** l = allocate2DArray(n,n);
    float** u = allocate2DArray(n,n);
    float*  x = (float*)malloc(n * sizeof(float));
    float*  y = (float*)malloc(n * sizeof(float));

    luDecomp(a, l, u, n);
    forwardSubstitution(l, b, y, n);
    backwardSubstitution(u, y, x, n);
    print1DMatrix(x,n);
}
//////////////////////////////////////////////////////////////////

int main(){
    int n = 3;
    float** a = allocate2DArray(n,n);
    float*  b = (float*)malloc(n * sizeof(float));
    float values[3][3] = {
        {1.0, 4.0, 5.0},
        {6.0, 8.0, 22.0},
        {32.0, 5.0, 5.0}
    };

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            a[i][j] = values[i][j];
        }
    }
    float bvalues[] = {1.0, 2.0, 3.0};
    for (int i = 0; i < n; i++) {
        b[i] = bvalues[i];
    }
    linalgsolve(a,b,n);
}