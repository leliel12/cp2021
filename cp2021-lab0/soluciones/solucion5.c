#include "wtime.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 500
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
    for (int y = 0; y < N; ++y) {
        /* NOTA: Lo unico que cambia es el orden de estos dos loops */
        for (int k = 0; k < N; ++k) {
            for (int x = 0; x < N; ++x) {
                c[x + y*N] += a[k + y*N] * b[x + k*N];
            }
        }
    }
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

