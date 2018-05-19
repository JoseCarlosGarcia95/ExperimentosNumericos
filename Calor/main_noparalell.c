#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define M_PI 3.14159265358979323846

float calcular_salto(float max, unsigned long n) {
    return max / n;
}

void resolver_ecuacion_calor_explicito(float xmax, float tmax, unsigned long n, unsigned long m,
                                       float c, float (*f)(float), float (*g)(float))
{
    unsigned long i, j;
    float dx, dt, **u, mu;

    dx = calcular_salto(xmax, n);
    dt = calcular_salto(tmax, m);

    u = malloc(sizeof(float*) * (n + 1));

    for (i = 0; i <= n; i++) {
        u[i] = malloc((m + 1) * sizeof(float));
    }

    for (j = 0; j <= m; j++) {
        u[0][j] = g(j * dt);
        u[n][j] = g(j * dt);
    }

    for (i = 0; i <= n; i++) {
        u[i][0] = f(i * dx);
    }

    mu = c * dt / (dx*dx);
    
    for (j = 0; j < m; j++) {
        
        for (i = 1; i < n; i++) {
            u[i][j+1] = u[i][j] + mu * (u[i-1][j] - 2*u[i][j] + u[i+1][j]);
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
    int n, m;

    if (argc > 2) {
        n = atoi(args[1]);
        m = atoi(args[2]);
        
        resolver_ecuacion_calor_explicito(4., 0.8, n, m, 0.005, f_ejemplo, g_ejemplo);
    }
}
