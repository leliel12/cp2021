#include "wtime.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#define N 500
#define STEPS 10

/* IMPORTANTE
 * ----------
 * La memoria es unidimensional, y en C las matrices multidimensionales se
 * representan de forma row-major: las filas se guardan de forma contigua.
 *
 * Matriz A 2x3:
 *  a11 a12 a13
 *  a21 a22 a23
 *
 * A en lenguajes row-major que cuentan desde 0 (C, C++):
 * Memoria: [a11 a12 a13 a21 a22 a23]
 * Indices:   0   1   2   3   4   5
 *
 * A en lenguajes column-major que cuentan desde 1 (Fortran, Julia):
 * Memoria: [a11 a21 a12 a22 a13 a23]
 * Indices:   1   2   3   4   5   6
 *
 */
void init_a(float * m)
{
    for (int y = 0; y < N; ++y) {
        for (int x = 0; x < N; ++x) {
            // m_{y,x} = y
            m[y * N + x] = y + 1.0f;
        }
    }
}


void init_b(float * m)
{
    for (int y = 0; y < N; ++y) {
        for (int x = 0; x < N; ++x) {
            // Triangulo superior U = 1, resto 0
            if (x < y) {
                m[y * N + x] = 0.0f;
            } else {
                m[y * N + x] = 1.0f;
            }
        }
    }
}

//inicializa b pero de manera transpuesta bt[x,y]=b[y,x]
void init_b_transposed(float * bt)
{
    for (int y = 0; y < N; ++y) {
        for (int x = 0; x < N; ++x) {
            // Triangulo superior U = 1, resto 0
            if (x <= y) {
                bt[y * N + x] = 1.0f;
            } else {
                bt[y * N + x] = 0.0f;
            }
        }
    }
}


void init_c(float * m)
{
    // C = 0
    for (int idx = 0; idx < N*N; ++idx) {
        m[idx] = 0.0f;
    }
}


bool check_result(const float * m)
{
    bool pass = true;
    for (int y = 0; y < N; ++y) {
        for (int x = 0; x < N; ++x) {
            float expected = (x+1) * (y+1) * STEPS;
            float got = m[y * N + x];
            if (got != expected) {
                printf("%d,%d: got %f, expected %f\n", y, x, got, expected);
                pass = false;
            }
        }
    }
    return pass;
}


void matmul_naive(const float * a, const float * b, float * c)
{
    /* FALTA: calcular C = A*B + C */
    for (int y = 0; y < N; ++y) {
        for (int x = 0; x < N; ++x) {
        	for (int k = 0; k < N; ++k) {
			c[y * N + x]+=a[y * N + k]* b[k * N + x];
		}
        }
    }

}

void matmul_naive_omp(const float * a, const float * b, float * c)
{
    /* FALTA: calcular C = A*B + C */
#pragma omp parallel for default(none) \
  	schedule(dynamic) \
  	shared(c,a,b)
    for (int y = 0; y < N; ++y) {
        for (int x = 0; x < N; ++x) {
        	for (int k = 0; k < N; ++k) {
			c[y * N + x]+=a[y * N + k]* b[k * N + x];
		}
        }
    }

}
// realiza el producto pero trasponiendo la matriz bt
// C = A*(B^T)^T + C
// lo que equivale a:
// C = A*B + C
void matmul_transposed(const float * a, const float * bt, float * c)
{
    /* FALTA: calcular C = A*B + C */
    for (int y = 0; y < N; ++y) {
        for (int x = 0; x < N; ++x) {
        	for (int k = 0; k < N; ++k) {
			// la idea es evitar que puntero sobre el que
			// itera k de saltos (evitar cache missing y TLB missing)
			c[y * N + x]+=a[y * N + k]* bt[x * N + k];
		}
        }
    }

}

void matmul_transposed_ptr(const float * a, const float * bt, float * c)
{
    /* FALTA: calcular C = A*B + C */
    for (int y = 0; y < N; ++y) {
        for (int x = 0; x < N; ++x) {
        	for (int k = 0; k < N; ++k) {
			// la idea es evitar que puntero sobre el que
			// itera k de saltos (evitar cache missing y TLB missing)
			*(c+ y*N+x)+=(*(a + y*N+ k))*(*(bt + x*N+ k));
		}
        }
    }

}


void matmul_transposed_omp(const float * a, const float * bt, float * c)
{
    /* FALTA: calcular C = A*B + C */
	  
#pragma omp parallel for default(none) \
  	schedule(dynamic) \
  	shared(c,a,bt)
    for (int y = 0; y < N; ++y) {
        for (int x = 0; x < N; ++x) {
        	for (int k = 0; k < N; ++k) {
			c[y * N + x]+=a[y * N + k]* bt[x * N + k];
		}
        }
    }

}

int main()
{

    size_t matsize = N * N * sizeof(float);

    /* pedir memoria para las matrices */
    float *a = malloc(matsize);
    float *b = malloc(matsize);
    float *c = malloc(matsize);


    /* inicializar valores */
    init_a(a);
    init_b(b);
    init_c(c);

    double start = wtime();
    for (int i = 0; i < STEPS; ++i) {
        matmul_naive(a, b, c);
    }

    double end = wtime();
    double elapsed = end - start;
    double operations = STEPS * (2.0 * N * N * N + N * N);
    double gflops = operations / (1000.0 * 1000.0 * 1000.0 * elapsed);
    if (check_result(c)) {
        printf("Naive run: %f GFLOPS\n", gflops);
    } else {
        printf("Resultado incorrecto!\n");
    }
    /* inicializar valores */
    init_a(a);
    init_b(b);
    init_c(c);

    start = wtime();
    for (int i = 0; i < STEPS; ++i) {
        matmul_naive_omp(a, b, c);
    }

    end = wtime();
    elapsed = end - start;
    gflops = operations / (1000.0 * 1000.0 * 1000.0 * elapsed);
    if (check_result(c)) {
        printf("Naive omp run: %f GFLOPS\n", gflops);
    } else {
        printf("Resultado incorrecto!\n");
    }



    /* inicializar valores */
    init_a(a);
    init_b_transposed(b);
    init_c(c);

    start = wtime();
    for (int i = 0; i < STEPS; ++i) {
        matmul_transposed(a, b, c);
    }

    end = wtime();
    elapsed = end - start;
    gflops = operations / (1000.0 * 1000.0 * 1000.0 * elapsed);
    if (check_result(c)) {
        printf("Transposed run: %f GFLOPS\n", gflops);
    } else {
        printf("Resultado incorrecto!\n");
    }

    /*init_a(a);
    init_b_transposed(b);
    init_c(c);

    start = wtime();
    for (int i = 0; i < STEPS; ++i) {
        matmul_transposed_ptr(a, b, c);
    }

    end = wtime();
    elapsed = end - start;
    gflops = operations / (1000.0 * 1000.0 * 1000.0 * elapsed);
    if (check_result(c)) {
        printf("Transposed Ptr run: %f GFLOPS\n", gflops);
    } else {
        printf("Resultado incorrecto!\n");
    }*/

    init_a(a);
    init_b_transposed(b);
    init_c(c);

    start = wtime();
    for (int i = 0; i < STEPS; ++i) {
        matmul_transposed_omp(a, b, c);
    }

    end = wtime();
    elapsed = end - start;
    gflops = operations / (1000.0 * 1000.0 * 1000.0 * elapsed);
    if (check_result(c)) {
        printf("Transposed OMP run: %f GFLOPS\n", gflops);
    } else {
        printf("Resultado incorrecto!\n");
    }


    /* devolver algun resultado para que el compilador no descarte codigo */
    return (int) c[0];
}

