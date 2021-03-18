#include "wtime.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>

#define N 3000
#define STEPS 10

void init(float * m) {
    for (int y = 0; y < N; ++y) {
        for (int x = 0; x < N; ++x) {
            /* Notar que las matrices son row-major */
            m[y * N + x] = 1.0f;
        }
    }
}

void matmul(const float * a, const float * b, float * c) {
    /*
     * Llamamos a una implementacion de BLAS (OpenBLAS, MKL, BLIS)
     * para comparar con el estado del arte.
     */
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                N, N, N,
                1.0, a, N,
                b, N,
                1.0, c, N);
}

int main(int argc, char ** argv) {

    size_t matsize = N * N * sizeof(float);

    /* alojar matrices */
    float *a = malloc(matsize);
    float *b = malloc(matsize);
    float *c = malloc(matsize);

    /* inicializar valores */
    init(a);
    init(b);
    init(c);

    double start = wtime();
    for (int i = 0; i < STEPS; ++i) {
        matmul(a, b, c);
    }

    double elapsed = wtime() - start;
    double operations = STEPS * (2.0 * N * N * N + N * N);
    double gflops = operations / (1000.0 * 1000.0 * 1000.0 * elapsed);
    printf("%f GFLOPS\n", gflops);

    /* devolver algun resultado para que el compilador no descarte codigo */
    return (int) c[0];
}

