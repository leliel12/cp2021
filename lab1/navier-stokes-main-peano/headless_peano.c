/*
  ======================================================================
   demo.c --- protoype to show off the simple solver
  ----------------------------------------------------------------------
   Author : Jos Stam (jstam@aw.sgi.com)
   Creation Date : Jan 9 2003

   Description:

	This code is a simple prototype that demonstrates how to use the
	code provided in my GDC2003 paper entitles "Real-Time Fluid Dynamics
	for Games". This code uses OpenGL and GLUT for graphics and interface

  =======================================================================
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "wtime.h"
#include "solver_peano.h"

/* global variables */

static int N;
static float dt, diff, visc;
static float force, source;

static struct float_link *u, *v, *u_prev, *v_prev;
static struct float_link *dens, *dens_prev;

/*
  ----------------------------------------------------------------------
   free/clear/allocate simulation data
  ----------------------------------------------------------------------
*/

static void free_data(void)
{
    if (u) {
        free(u);
    }
    if (v) {
        free(v);
    }
    if (u_prev) {
        free(u_prev);
    }
    if (v_prev) {
        free(v_prev);
    }
    if (dens) {
        free(dens);
    }
    if (dens_prev) {
        free(dens_prev);
    }
}

static void clear_data(void)
{
    int i, size = (N + 2) * (N + 2);

    for (i = 0; i < size; i++) {
        u[i].val[0] = v[i].val[0] = u_prev[i].val[0] = v_prev[i].val[0] = dens[i].val[0] = dens_prev[i].val[0] = 0.0f;
        u[i].val[1] = v[i].val[1] = u_prev[i].val[1] = v_prev[i].val[1] = dens[i].val[1] = dens_prev[i].val[1] = 0.0f;
    }
}

static int allocate_data()
{
    int size   = (N + 2) * (N + 2);
    int sizein = N*N;
    u         = (struct float_link *)malloc(size * sizeof(struct float_link));
    v         = (struct float_link *)malloc(size * sizeof(struct float_link));
    u_prev    = (struct float_link *)malloc(size * sizeof(struct float_link));
    v_prev    = (struct float_link *)malloc(size * sizeof(struct float_link));
    dens      = (struct float_link *)malloc(size * sizeof(struct float_link));
    dens_prev = (struct float_link *)malloc(size * sizeof(struct float_link));
    
    init_pointers(u        ,N);
    init_pointers(v        ,N);
    init_pointers(u_prev   ,N);
    init_pointers(v_prev   ,N);
    init_pointers(dens     ,N);
    init_pointers(dens_prev,N);
	      
    if (!u || !v || !u_prev || !v_prev || !dens || !dens_prev) {
        fprintf(stderr, "cannot allocate data\n");
        return (0);
    }

    //fclose(fp);
    //exit(-1);
    return (1);
}

static void react(struct float_link *d,struct float_link *u,struct float_link *v)
{
    int i ;
    int size = (N + 2) * (N + 2);
    int sizen = N*N;
    float max_velocity2 = 0.0f;
    float max_density = 0.0f;

    max_velocity2 = max_density = 0.0f;
    for (i = 0; i < size; i++) {
        if (max_velocity2 < u[i].val[0] * u[i].val[0] + v[i].val[0] * v[i].val[0]) {
            max_velocity2 = u[i].val[0] * u[i].val[0] + v[i].val[0] * v[i].val[0];
        }
        if (max_density < d[i].val[0]) {
            max_density = d[i].val[0];
        }
    }

    for (i = 0; i < size; i++) {
        u[i].val[0] = v[i].val[0] = d[i].val[0] = 0.0f;
    }

    if (max_velocity2 < 0.0000005f) {
        u[sizen/2].val[0] = force * 10.0f;
        v[sizen/2].val[0] = force * 10.0f;
    }
    if (max_density < 1.0f) {
        d[sizen/2].val[0] = source * 10.0f;
    }

    return;
}

static void one_step(int step, int write_data)
{
    static int times = 1;
    static double start_t = 0.0;
    static double one_second = 0.0;
    static double react_ns_p_cell = 0.0;
    static double vel_ns_p_cell = 0.0;
    static double dens_ns_p_cell = 0.0;
    struct float_link *x, *x0;

    start_t = wtime();
    react(dens_prev, u_prev, v_prev);
    react_ns_p_cell += 1.0e9 * (wtime() - start_t) / (N * N);
    
    start_t = wtime();
    vel_step(N, u, v, u_prev, v_prev, visc, dt);
    vel_ns_p_cell += 1.0e9 * (wtime() - start_t) / (N * N);
    
    start_t = wtime();
    dens_step(N, dens, dens_prev, u, v, diff, dt);
    dens_ns_p_cell += 1.0e9 * (wtime() - start_t) / (N * N);

    if(write_data)
    {
        char buffer[10]="";
        char filename[50]="";
        strcpy(filename,"results/results-headless-peano-");
        sprintf(buffer, "%04d", step);
        strcat(filename,buffer);
        strcat(filename,".dat");
        printf("%s\n",filename);
        FILE *fp;
        fp = fopen(filename, "w");
        for (int i = 0; i <(N+2)*(N+2); i++)
        {
                int ii,jj;
                peano_keyto2d(N+2,i,&ii, &jj); 
                fprintf(fp, "%d %d %e %e %e\n",ii,jj,dens[i].val[0],u[i].val[0],v[i].val[0]);
        }
        fclose(fp);
    }
    else
    {
        if (1.0 < wtime() - one_second) { /* at least 1s between stats */
            printf("%lf, %lf, %lf, %lf: ns per cell total, react, vel_step, dens_step\n",
                   (react_ns_p_cell + vel_ns_p_cell + dens_ns_p_cell) / times,
                   react_ns_p_cell / times, vel_ns_p_cell / times, dens_ns_p_cell / times);
            fflush(stdout);
            one_second = wtime();
            react_ns_p_cell = 0.0;
            vel_ns_p_cell = 0.0;
            dens_ns_p_cell = 0.0;
            times = 1;
        } else {
            times++;
        }
    }



}


/*
  ----------------------------------------------------------------------
   main --- main routine
  ----------------------------------------------------------------------
*/

int main(int argc, char** argv)
{
    int i = 0;
    int b = 0;

    if (argc != 1 && argc != 7 && argc !=2) {
        fprintf(stderr, "usage : %s N dt diff visc force source\n", argv[0]);
        fprintf(stderr, "where:\n");
        fprintf(stderr, "\t exponent b : N=2^b, where N is grid resolution\n");
        fprintf(stderr, "\t dt     : time step\n");
        fprintf(stderr, "\t diff   : diffusion rate of the density\n");
        fprintf(stderr, "\t visc   : viscosity of the fluid\n");
        fprintf(stderr, "\t force  : scales the mouse movement that generate a force\n");
        fprintf(stderr, "\t source : amount of density that will be deposited\n");
        exit(1);
    }

    if (argc == 1) {
        N = 126; b=7;
        dt = 0.1f;
        diff = 0.0f;
        visc = 0.0f;
        force = 5.0f;
        source = 100.0f;
        fprintf(stderr, "Using defaults : N=%d dt=%g diff=%g visc=%g force = %g source=%g\n",
                N, dt, diff, visc, force, source);
    } 
    if (argc == 2) {
        b = atoi(argv[1]);
	N = 1 << b;
	N = N -2;
        dt = 0.1f;
        diff = 0.0f;
        visc = 0.0f;
        force = 5.0f;
        source = 100.0f;
        fprintf(stderr, "Using defaults : N=%d dt=%g diff=%g visc=%g force = %g source=%g\n",
                N, dt, diff, visc, force, source);
    }else {
        b = atoi(argv[1]);
	N = 1 << b;
	N= N-2;
        dt = atof(argv[2]);
        diff = atof(argv[3]);
        visc = atof(argv[4]);
        force = atof(argv[5]);
        source = atof(argv[6]);
    }

    if (!allocate_data()) {
        exit(1);
    }
    clear_data();

    for (i = 0; i < 2048; i++) {
        one_step(i,1);
    }

   free_data();
    exit(0);
}
