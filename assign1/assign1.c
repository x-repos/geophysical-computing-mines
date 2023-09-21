#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


int main() {
    int nt = 121;
    float dt = 0.02;

    float xmin = 0.0, xmax = 1.0, ymin = 0.0, ymax = 1.0;
    int nx = 101, ny = 101;
    float dx = (xmax - xmin) / (nx - 1);
    float dy = (ymax - ymin) / (ny - 1);
    float sigma = 0.05;
    float mx = (xmax - xmin) / 2;
    float my = (ymax - ymin) / 2;
    float vx = -0.1;
    float vy = 0.1;
    float L = 0.0;
    float** U = (float**)malloc(nx * sizeof(float*));
    float** Utmp = (float**)malloc(nx * sizeof(float*));
    for (int i = 0; i < nx; i++) {
        U[i] = (float*)malloc(ny * sizeof(float));
        for (int j = 0; j < ny; j++) {
            float xxx = xmin + i * dx;
            float yyy = ymin + j * dy;
            U[i][j] = exp(-(pow(xxx - mx, 2) / (2 * pow(sigma, 2)) + pow(yyy - my, 2) / (2 * pow(sigma, 2))));
        }
    }

    void Advection_Lax_Wendroff_Solution(float **u, float vx, float vy, float dx, float dy, float dt, float L, int nx, int ny) {
        float **u_temp = (float**)malloc(nx * sizeof(float*));

        // Allocate memory for temporary u
        for (int i = 0; i < nx; i++) {
            u_temp[i] = (float*)malloc(ny * sizeof(float));
        }

        float A, B;
        float term1, term2, term3, term4, term5;

        for (int i = 1; i < nx - 1; i++) {
            for (int j = 1; j < ny - 1; j++) {
                A = -vy;
                B = vx;
                term1 = -dt / (2 * dx) * A * (u[i + 1][j] - u[i - 1][j]);
                term2 = -dt / (2 * dy) * B * (u[i][j + 1] - u[i][j - 1]);
                term3 = (dt * dt) / (2 * dx * dx) * A * A * (u[i - 1][j] - 2 * u[i][j] + u[i + 1][j]);
                term4 = (dt * dt) / (2 * dy * dy) * B * B * (u[i][j - 1] - 2 * u[i][j] + u[i][j + 1]);
                term5 = (dt * dt) / (8 * dx * dy) * (A * B + B * A) * (
                    (u[i + 1][j + 1] - u[i - 1][j + 1]) - (u[i + 1][j - 1] - u[i - 1][j - 1])
                );

                u_temp[i][j] = u[i][j] + term1 + term2 + term3 + term4 + term5;
            }
        }

        for (int i = 1; i < nx - 1; i++) {
            for (int j = 1; j < ny - 1; j++) {
                u[i][j] = u_temp[i][j];
            }
        }

        // Boundary conditions
        for (int i = 0; i < nx; i++) {
            u[i][0] = L;
            u[i][ny - 1] = L;
        }
        for (int j = 0; j < ny; j++) {
            u[0][j] = L;
            u[nx - 1][j] = L;
        }
        for (int i = 0; i < nx; i++) {
            free(u_temp[i]);
        }
        free(u_temp);
    }


    clock_t start_time = clock();
    
    for (int i = 0; i < nt; i++)
    {   
        Advection_Lax_Wendroff_Solution(U, vx, vy, dx, dy, dt, L, nx, ny);


    }
    clock_t end_time = clock();
    double elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("C code execution time: %lf s\n", elapsed_time);

    FILE *file = fopen("U_matrix.txt", "w");

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            fprintf(file, "%lf ", U[i][j]);
        }
        fprintf(file, "\n"); 
    }
    fclose(file);

    // Free memory
    for (int i = 0; i < nx; i++) {
        free(U[i]);
    }
    free(U);
    return 0;
}