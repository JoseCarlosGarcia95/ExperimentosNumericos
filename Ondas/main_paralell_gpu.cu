#include <stdio.h>
#include <stdlib.h>
#include <math.h>

float calcular_salto(float max, unsigned long n) {
    return max / n;
}

__global__ void iteracion_gpu(float * u, unsigned long j, unsigned long m, float lambda) {
    unsigned int i;

    i = threadIdx.x + blockIdx.x * blockDim.x + 1;

    u[i*m + j + 1] = 2*(1-lambda*lambda)*u[i*m + j] + lambda * lambda  *(u[(i+1)*m + j] + u[(i-1)*m + j]) - u[i*m + j-1];
}

void resolver_ecuacion_onda_explicito(float xmax, float tmax, unsigned long n, unsigned long m,
                                       float c, float (*f)(float), float (*g)(float))
{
    unsigned long i, j;
    float dx, dt, *u, lambda, *u_gpu;
    cudaError_t cuda_status;

    dx = calcular_salto(xmax, n);
    dt = calcular_salto(tmax, m);

    u = (float*)malloc(sizeof(float) * (n + 1) * (m + 1));
    cuda_status = cudaMalloc(&u_gpu, sizeof(float) * (n + 1) * (m + 1));
    
    if (cuda_status != cudaSuccess) {
        printf("cudaMalloc returned error code %d\n", cuda_status);
        return;
    }
    
    for (j = 1; j <= m; j++) {
        u[j] = 0;
        u[n*m + j] = 0;
    }

    u[0] = f(0);
    u[n*m] = f(xmax);

    lambda = c * dt / dx;

    for (i = 1; i < n; i++) {
        u[i*m] = f(i * dx);
        u[i*m + 1] = (1-lambda*lambda)*f(i * dx) + lambda*lambda/2 * (f((i+1)*dx) + f((i-1)*dx)) + dt * g(i * dx);
    }
    
    cudaMemcpy(u_gpu, u, (n + 1) * (m + 1) * sizeof(float), cudaMemcpyHostToDevice);

    for (j = 1; j < m; j++) {        
        iteracion_gpu<<<1, n-1>>>(u_gpu, j, m, lambda);

        cuda_status = cudaDeviceSynchronize();
		if (cuda_status != cudaSuccess) {
			printf("cudaDeviceSynchronize returned error code %d\n", cuda_status);
            return;
        }
    }
    
    cudaMemcpy(u, u_gpu, (n + 1) * (m + 1) * sizeof(float), cudaMemcpyDeviceToHost);

    cudaFree(u_gpu);    
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

        resolver_ecuacion_onda_explicito(4., 0.8, n, m, 0.005, f_ejemplo, g_ejemplo);
    }
}
