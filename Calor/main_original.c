#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <pthread.h>

#define M_PI 3.14159265358979323846

double calcular_salto(double max, unsigned long n) {
    return max / n;
}

void dibujar_grafica(unsigned long puntos, double **grafica) {
    int i;
    FILE* gnuplot_pipe;
    gnuplot_pipe = popen ("gnuplot -persistent", "w");
    fprintf(gnuplot_pipe, "set title 'Gráfica de la función'\n");
    fprintf(gnuplot_pipe, "splot '-'\n");
    for (i = 0; i < puntos; i++) {
        fprintf(gnuplot_pipe, "%f %f %f\n", grafica[i][0], grafica[i][1], grafica[i][2]);
    }
    fprintf(gnuplot_pipe, "e\n");
    fprintf(gnuplot_pipe, "refresh\n");
    fclose(gnuplot_pipe);
}

void resolver_ecuacion_calor_explicito(double xmax, double tmax, unsigned long n, unsigned long m,
                                       double c, double (*f)(double), double (*g)(double))
{
    unsigned long k, i, j;
    double dx, dt, **u, mu, **grafica;

    dx = calcular_salto(xmax, n);
    dt = calcular_salto(tmax, m);

    u = malloc(sizeof(double*) * (n + 1));

    for (i = 0; i <= n; i++) {
        u[i] = malloc((m + 1) * sizeof(double));
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
    
    
    k = 0;
    /*
    grafica = malloc(sizeof(double*) * (n + 1) * (m + 1));

    for (i = 0; i <= n; i++) {
        for (int j = 0; j <= m; j++) {
            grafica[k] = malloc(sizeof(double) * 3);
            grafica[k][0]   = i * dx;
            grafica[k][1]   = j * dt;
            grafica[k][2]   = u[i][j];

            k++;
        }
    }
    
    dibujar_grafica((n+1)*(m+1), grafica);

    k = 0;
    for (i = 0; i <= n; i++) {
        free(u[i]);
        for (int j = 0; j <= m; j++) {
            free(grafica[k++]);
        }
    }

    free(grafica);*/
    free(u);
}

double f_ejemplo(double x) {
    return sin(M_PI * x / 4) * (1 + 2 * cos(M_PI * x / 4));
}

double g_ejemplo(double t) {
    return 0;
}

int main(int argc, char** args) {
    resolver_ecuacion_calor_explicito(4., 0.8, 20, 20, 0.005, f_ejemplo, g_ejemplo);
}
