#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

#define M_PI 3.14159265358979323846

float calcular_salto(float max, unsigned long n) {
    return max / n;
}

void resolver_ecuacion_onda_explicito(float xmax, float tmax, unsigned long n, unsigned long m,
                                       float c, float (*f)(float), float (*g)(float))
{
    unsigned long i, j;
    float dx, dt, **u, lambda;

    dx = calcular_salto(xmax, n);
    dt = calcular_salto(tmax, m);

    u = malloc(sizeof(float*) * (n + 1));

    for (i = 0; i <= n; i++) {
        u[i] = malloc((m + 1) * sizeof(float));
    }

    for (j = 1; j <= m; j++) {
        u[0][j] = 0;
        u[n][j] = 0;
    }

    u[0][0] = f(0);
    u[n][0] = f(xmax);
    
    lambda = c * dt / dx;

    for (i = 1; i < n; i++) {
        u[i][0] = f(i * dx);
        u[i][1] = (1-lambda*lambda)*f(i * dx) + lambda*lambda/2 * (f((i+1)*dx) + f((i-1)*dx)) + dt * g(i * dx);
    }
    
    for (j = 1; j < m; j++) {
        
#pragma omp parallel for
        for (i = 1; i < n; i++) {
            u[i][j+1] = 2*(1-lambda*lambda)*u[i][j] + lambda*lambda*(u[i+1][j] + u[i-1][j]) - u[i][j-1];
        }
        
    }

    for (i = 0; i <= n; i++) {
        free(u[i]);
    }
    
    free(u);
}

float f_ejemplo(float x) {
    return sin(M_PI * x / 4) * (1 + 2 * cos(M_PI * x / 4));
}

float g_ejemplo(float t) {
    return 0;
}

int main(int argc, char** args) {
    int n, m, hilos;

    if (argc > 3) {
        n = atoi(args[1]);
        m = atoi(args[2]);
        hilos = atoi(args[3]);

        omp_set_num_threads(hilos);

        resolver_ecuacion_onda_explicito(4., 0.8, n, m, 0.005, f_ejemplo, g_ejemplo);
    }
    
    return 0; 
}
