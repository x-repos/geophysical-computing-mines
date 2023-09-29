


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// Allocation
//////////////////////////////////////////////////////////////////
float* allocate1DArray(int n){
    /*
    float* a = allocate1DArray(n);
    */
    float* array = (float*)malloc(n * sizeof(float));
    for (int i = 0; i < n; i++)
    {
        array[i] = 0.0;
    }
    
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


void free2DArray(float** array, int rows) {
    for (int i = 0; i < rows; i++) {
        free(array[i]);
    }
    free(array);
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

float* linalgsolve(float** a, float* b, int n){
    float** l = allocate2DArray(n,n);
    float** u = allocate2DArray(n,n);
    float*  x = (float*)malloc(n * sizeof(float));
    float*  y = (float*)malloc(n * sizeof(float));

    luDecomp(a, l, u, n);
    forwardSubstitution(l, b, y, n);
    backwardSubstitution(u, y, x, n);
    free2DArray(l,n);
    free2DArray(u,n);
    return x;
}


// Linspace

void linspace(float start, float end, int num, float* result) {
    float delta = (end - start) / (num - 1);
    for (int i = 0; i < num; i++) {
        result[i] = start + delta * i;
    }
}

// save matrix
//////////////////////////////////////////////////////////////////
void saveMatrixToFile(float **matrix, int rows, int cols, const char *filename, const char *mode) {
    FILE *file = fopen(filename, mode);
    if (file == NULL) {
        printf("Error opening file!\n");
        return;
    }

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            fprintf(file, "%lf ", matrix[i][j]); // write each element with a space in between
        }
        fprintf(file, "\n"); // new line for each row
    }

    fclose(file);
}

void saveArrayToFile(float *data, int n, const char *filename) {
    FILE *fp = fopen(filename, "w"); // Open the file for writing
    if (fp == NULL) {
        printf("Error opening file\n");
        return; // Exit if file can't be opened
    }

    for (int i = 0; i < n; i++) {
        fprintf(fp, "%lf\n", data[i]); // Write each element to the file, each on a new line
    }

    fclose(fp); // Close the file
}


// Internal function
void Setup_Tridiagonal_Heterogeneous(float** A, float* a, float* b, float* c, int n) {
    int i;
    for (i = 0; i < n; i++) {
        A[i][i] = b[i];
    }
    for (i = 1; i < n; i++) {
        A[i][i - 1] = a[i - 1];
    }
    for (i = 0; i < n - 1; i++) {
        A[i][i + 1] = c[i];
    }
    A[0][0] = 1.0;
    A[0][1] = 0.0;
    A[n - 1][n - 1] = 1.0;
    A[n - 1][n - 2] = 0.0;
}

void ADI_Solution_Heterogeneous(float **U, float** K, float** F, float h, int n){
    int i, j, ii, jj; 
    int nx = n, ny = n;
    float** U1 = allocate2DArray(n,n);
    for (i = 0; i < nx; i++) {
        U1[i][0] = U[i][0];
    }
    for (j = 0; j < ny; j++) {
        U1[0][j] = U[0][j];
    }
    for (i = 0; i < nx; i++) {
        U1[i][ny-1] = U[i][ny-1];
    }
    for (j = 0; j < ny; j++) {
        U1[nx-1][j] = U[nx-1][j];
    }

    float *axmat, *bxmat, *cxmat, *bx, *ktemp;
    float *aymat, *bymat, *cymat, *by; 
    axmat = allocate1DArray(nx-1);
    bxmat = allocate1DArray(nx);
    cxmat = allocate1DArray(nx-1);
    bx    = allocate1DArray(nx);
    aymat = allocate1DArray(ny-1);
    bymat = allocate1DArray(ny);
    cymat = allocate1DArray(ny-1);
    by    = allocate1DArray(ny);
    ktemp = allocate1DArray(nx);

    for (jj = 1; jj < ny-1; jj++) {
        for (i = 0; i < nx-1; i++) {
            axmat[i] = K[i][jj];
        }
        for (int i = 1; i < nx; i++)
        {
            ktemp[i] = K[i-1][jj];
        }
        for (i = 0; i < nx; i++) {
            bxmat[i] = -2.0 * K[i][jj] - K[i][jj-1]-ktemp[i];
        }
        for (i = 1; i < nx; i++) {
            cxmat[i-1] = K[i][jj];
        }
        float** Ax = allocate2DArray(n,n);
        Setup_Tridiagonal_Heterogeneous(Ax, axmat, bxmat, cxmat, nx);

        for (i = 0; i < nx; i++) {
            bx[i] = -h * h * F[i][jj] - K[i][jj] * U[i][jj+1] - K[i][jj-1] * U[i][jj-1];
        }
        bx[0] = U[0][jj];
        bx[nx-1] = U[nx-1][jj];

        float *result = linalgsolve(Ax, bx, nx);
        for (i = 0; i < nx; i++) {
            U1[i][jj] = result[i];
        }
        free(result);
        for (int i = 0; i < nx; i++) {
            free(Ax[i]);
        }
        free(Ax);        
    }

    for (ii = 1; ii < nx-1; ii++) {

        for (j = 0; j < ny-1; j++) {
            aymat[j] = K[ii][j];
        }
        for (int j = 1; j < ny; j++)
        {
            ktemp[j] = K[ii][j-1];
        }
        for (j = 0; j < ny; j++) {
            bymat[j] = -2.0 * K[ii][j] - K[ii-1][j] - ktemp[j];
        }
        for (j = 1; j < ny; j++) {
            cymat[j-1] = K[ii][j];
        }
        float** Ay = allocate2DArray(n,n);
        Setup_Tridiagonal_Heterogeneous(Ay, aymat, bymat, cymat, ny);
        for (j = 0; j < ny; j++) {
            by[j] = -h * h * F[ii][j] - K[ii][j] * U1[ii+1][j] - K[ii-1][j] * U1[ii-1][j];
        }
        by[0] = U1[ii][0];
        by[ny-1] = U1[ii][ny-1];
        float *result = linalgsolve(Ay, by, ny);
        for (j = 0; j < ny; j++) {
            U[ii][j] = result[j];
        }

        free(result);
        for (int i = 0; i < nx; i++) {
            free(Ay[i]);
        }
        free(Ay);
       
    }
    for (int i = 0; i < nx; i++) {
        free(U1[i]);
    }
    free(U1);
}


// Main
//////////////////////////////////////////////////////////////////
int main(){
    float t1 = -1.0, t2 = 1.0, a = 1, b = 1;
    float xmin = 0.0, ymin = 0.0, xmax, ymax;
    xmax = a, ymax = b;
    int n = 100, nt = 500;
    float* xx = allocate1DArray(n);
    float* yy = allocate1DArray(n);
    linspace(xmin, xmax, n, xx);
    linspace(ymin, ymax, n, yy);
    float h = xx[1] - xx[0];
    float** u0 = allocate2DArray(n,n);

    for (int j = 0; j < n; j++) {
        u0[0][j] = t1;
    }
    for (int j = 0; j < n; j++) {
        u0[n-1][j] = 0.0;
    }
    for (int i = 0; i < n; i++) {
        u0[i][0] = 0.0;
    }
    for (int i = 0; i < n; i++) {
        u0[i][n-1] = t2;
    }
    u0[0][0] = 0.0;
    u0[0][n-1] = 0.0;
    u0[n-1][0] = 0.0;
    u0[n-1][n-1] = 0.0;

    float** K = allocate2DArray(n,n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            K[i][j] = 1.0 - 0.75 * exp(- (xx[j] - a/2.0) * (xx[j] - a/2.0) - (yy[i] - b/2.0) * (yy[i] - b/2.0));
        }
    }
    float sourcesize = 0.1;
    int halfsize = (int) round(sourcesize / a / h / 2.0);

    float** ff = allocate2DArray(n,n);

    for (int i = (n/2) - 1 - halfsize; i < (n/2) - 1 + halfsize; i++) {
        for (int j = (n/2) - 1 - halfsize; j < (n/2) - 1 + halfsize; j++) {
            ff[i][j] = 200.0 * t2;
        }
    }

    for (int i = 0; i < nt; i++)
    {
        ADI_Solution_Heterogeneous(u0,K,ff,h,n);
        saveMatrixToFile(u0,n,n,"results.txt","a");
    }
    free2DArray(K,n);
    free2DArray(ff,n);
}
