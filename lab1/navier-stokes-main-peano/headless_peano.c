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
        u[i].val = v[i].val = u_prev[i].val = v_prev[i].val = dens[i].val = dens_prev[i].val = 0.0f;
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
    //FILE *fp;
    //fp = fopen("mem_order.dat", "w");//opening file.
    //inside
/*
    for (i = 0; i <N; i++) 
    {
       for (j = 0; j <N; j++) 
       {
	       k  = peano_2dtokey(N,i,j);
	       if(j==(N-1))
	       {
		   kU = sizein + i;//borde superior
	       }
	       else 
	       {
	           kU = peano_2dtokey(N,i,j+1);
	       }
	       if(j==0)
	       {
		   kD = sizein + N + i;//borde inferior
	       }
	       else
	       {
		   kD = peano_2dtokey(N,i,j-1);
	       }
	       if(i==(N-1))
	       {
	           kR = sizein + 2*N + j  ;//borde derecho
	       }
	       else
	       {
	           kR = peano_2dtokey(N,i+1,j);
	       }
	       if(i==0)
	       {
		   kL =  sizein + 3*N + j;//borde izquierdo
	       }
	       else
	       {
	           kL = peano_2dtokey(N, i-1, j   );
	       }

               //fprintf(fp, "%d %d %d %d %d %d %d\n",i,j,k,kU,kD,kL,kR);
	       u[k].x = i;
	       u[k].y = j;
	       u[k].U = &u[kU];
	       u[k].D = &u[kD];
	       u[k].L = &u[kL];
	       u[k].R = &u[kR];

	       v[k].x = i;
	       v[k].y = j;
	       v[k].U = &v[kU];
	       v[k].D = &v[kD];
	       v[k].L = &v[kL];
	       v[k].R = &v[kR];

	       u_prev[k].x = i;
	       u_prev[k].y = j;
	       u_prev[k].U = &u_prev[kU];
	       u_prev[k].D = &u_prev[kD];
	       u_prev[k].L = &u_prev[kL];
	       u_prev[k].R = &u_prev[kR];

	       v_prev[k].x = i;
	       v_prev[k].y = j;
	       v_prev[k].U = &v_prev[kU];
	       v_prev[k].D = &v_prev[kD];
	       v_prev[k].L = &v_prev[kL];
	       v_prev[k].R = &v_prev[kR];

	       dens[k].x = i;
	       dens[k].y = j;
	       dens[k].U = &dens[kU];
	       dens[k].D = &dens[kD];
	       dens[k].L = &dens[kL];
	       dens[k].R = &dens[kR];

	       dens_prev[k].x = i;
	       dens_prev[k].y = j;
	       dens_prev[k].U = &dens_prev[kU];
	       dens_prev[k].D = &dens_prev[kD];
	       dens_prev[k].L = &dens_prev[kL];
	       dens_prev[k].R = &dens_prev[kR];

       }
    }

    //boundaries
    i=0;
    for (int t = 0; t <N; t++) 
    {
	       // borde superior j=N
	       s = sizein + t;
	       u[s].x = t;
	       u[s].y = N;
	       v[s].x = u_prev[s].x = v_prev[s].x= v_prev[s].x=dens[s].x= dens_prev[s].x =u[s].x;
	       v[s].y = u_prev[s].y = v_prev[s].y= v_prev[s].y=dens[s].y= dens_prev[s].y =u[s].y;

	       //apunta para abajo
	       k  = peano_2dtokey(N,t,N-1);
               //fprintf(fp, "%d %d %d %d %d %d %d\n",t,N,s,k,k,k,k);
	       u[s].U         = u[s].D         = u[s].L         = u[s].R         = &u[k];
	       v[s].U         = v[s].D         = v[s].L         = v[s].R         = &v[k];
	       u_prev[s].U    = u_prev[s].D    = u_prev[s].L    = u_prev[s].R    = &u_prev[k];
	       v_prev[s].U    = v_prev[s].D    = v_prev[s].L    = v_prev[s].R    = &v_prev[k];
	       dens[s].U      = dens[s].D      = dens[s].L      = dens[s].R      = &dens[k];
	       dens_prev[s].U = dens_prev[s].D = dens_prev[s].L = dens_prev[s].R = &dens_prev[k];
	       // borde inferior j=-1
	       s = sizein + N+ t;
	       u[s].x =  t;
	       u[s].y = -1;
	       v[s].x = u_prev[s].x = v_prev[s].x= v_prev[s].x=dens[s].x= dens_prev[s].x =u[s].x;
	       v[s].y = u_prev[s].y = v_prev[s].y= v_prev[s].y=dens[s].y= dens_prev[s].y =u[s].y;
	       //apunta para arriba
	       k  = peano_2dtokey(N,t,0);
               //fprintf(fp, "%d %d %d %d %d %d %d\n",t,-1,s,k,k,k,k);
	       u[s].U         = u[s].D         = u[s].L         = u[s].R         = &u[k];
	       v[s].U         = v[s].D         = v[s].L         = v[s].R         = &v[k];
	       u_prev[s].U    = u_prev[s].D    = u_prev[s].L    = u_prev[s].R    = &u_prev[k];
	       v_prev[s].U    = v_prev[s].D    = v_prev[s].L    = v_prev[s].R    = &v_prev[k];
	       dens[s].U      = dens[s].D      = dens[s].L      = dens[s].R      = &dens[k];
	       dens_prev[s].U = dens_prev[s].D = dens_prev[s].L = dens_prev[s].R = &dens_prev[k];
   	       // borde derecho i=N
	       s = sizein + 2*N+ t;
	       u[s].x = N;
	       u[s].y = t;
	       v[s].x = u_prev[s].x = v_prev[s].x= v_prev[s].x=dens[s].x= dens_prev[s].x =u[s].x;
	       v[s].y = u_prev[s].y = v_prev[s].y= v_prev[s].y=dens[s].y= dens_prev[s].y =u[s].y;
	       //apunta a izquierda
	       k  = peano_2dtokey(N,N-1,t);
               //fprintf(fp, "%d %d %d %d %d %d %d\n",N,t,s,k,k,k,k);
	       u[s].U         = u[s].D         = u[s].L         = u[s].R         = &u[k];
	       v[s].U         = v[s].D         = v[s].L         = v[s].R         = &v[k];
	       u_prev[s].U    = u_prev[s].D    = u_prev[s].L    = u_prev[s].R    = &u_prev[k];
	       v_prev[s].U    = v_prev[s].D    = v_prev[s].L    = v_prev[s].R    = &v_prev[k];
	       dens[s].U      = dens[s].D      = dens[s].L      = dens[s].R      = &dens[k];
	       dens_prev[s].U = dens_prev[s].D = dens_prev[s].L = dens_prev[s].R = &dens_prev[k];
   	       // borde izquierdo i=-1
	       s = sizein + 3*N+ t;
	       u[s].x = -1;
	       u[s].y =  t;
	       v[s].x = u_prev[s].x = v_prev[s].x= v_prev[s].x=dens[s].x= dens_prev[s].x =u[s].x;
	       v[s].y = u_prev[s].y = v_prev[s].y= v_prev[s].y=dens[s].y= dens_prev[s].y =u[s].y;
	       //apunta a derecha
	       k  = peano_2dtokey(N,0,t);
               //fprintf(fp, "%d %d %d %d %d %d %d\n",-1,t,s,k,k,k,k);
	       u[s].U         = u[s].D         = u[s].L         = u[s].R         = &u[k];
	       v[s].U         = v[s].D         = v[s].L         = v[s].R         = &v[k];
	       u_prev[s].U    = u_prev[s].D    = u_prev[s].L    = u_prev[s].R    = &u_prev[k];
	       v_prev[s].U    = v_prev[s].D    = v_prev[s].L    = v_prev[s].R    = &v_prev[k];
	       dens[s].U      = dens[s].D      = dens[s].L      = dens[s].R      = &dens[k];
	       dens_prev[s].U = dens_prev[s].D = dens_prev[s].L = dens_prev[s].R = &dens_prev[k];
    }
    // esquinas
       // esquina abajo izquierda
       s = sizein + 4*N;
       u[s].x = -1;
       u[s].y = -1;
       v[s].x = u_prev[s].x = v_prev[s].x= v_prev[s].x=dens[s].x= dens_prev[s].x =u[s].x;
       v[s].y = u_prev[s].y = v_prev[s].y= v_prev[s].y=dens[s].y= dens_prev[s].y =u[s].y;
       kL = sizein + N;   //sobre el borde a la derecha de la esquina
       kU = sizein + 3*N; //sobre el borde arriba de la esquina
       u[s].U         = u[s].D         = &u[kU];        
       u[s].L         = u[s].R         = &u[kL];
       v[s].U         = v[s].D         = &v[kU];
       v[s].L         = v[s].R         = &v[kL];
       u_prev[s].U    = u_prev[s].D    = &u_prev[kU];
       u_prev[s].L    = u_prev[s].R    = &u_prev[kL];
       v_prev[s].U    = v_prev[s].D    = &v_prev[kU];
       v_prev[s].L    = v_prev[s].R    = &v_prev[kL];
       dens[s].U      = dens[s].D      = &dens[kU];
       dens[s].L      = dens[s].R      = &dens[kL];
       dens_prev[s].U = dens_prev[s].D = &dens_prev[kU];
       dens_prev[s].L = dens_prev[s].R = &dens_prev[kL];

      // esquina abajo derecha
       s = sizein + 4*N +1;
       u[s].x = N;
       u[s].y = -1;
       v[s].x = u_prev[s].x = v_prev[s].x= v_prev[s].x=dens[s].x= dens_prev[s].x =u[s].x;
       v[s].y = u_prev[s].y = v_prev[s].y= v_prev[s].y=dens[s].y= dens_prev[s].y =u[s].y;
       kL = sizein + 2*N-1; //sobre el borde a la izquierda de la esquina
       kU = sizein + 2*N;   //sobre el borde arriba de la esquina
       u[s].U         = u[s].D         = &u[kU];        
       u[s].L         = u[s].R         = &u[kL];
       v[s].U         = v[s].D         = &v[kU];
       v[s].L         = v[s].R         = &v[kL];
       u_prev[s].U    = u_prev[s].D    = &u_prev[kU];
       u_prev[s].L    = u_prev[s].R    = &u_prev[kL];
       v_prev[s].U    = v_prev[s].D    = &v_prev[kU];
       v_prev[s].L    = v_prev[s].R    = &v_prev[kL];
       dens[s].U      = dens[s].D      = &dens[kU];
       dens[s].L      = dens[s].R      = &dens[kL];
       dens_prev[s].U = dens_prev[s].D = &dens_prev[kU];
       dens_prev[s].L = dens_prev[s].R = &dens_prev[kL];


       // esquina arriba izquierda
       s = sizein + 4*N + 2;
       u[s].x = -1;
       u[s].y = N;
       v[s].x = u_prev[s].x = v_prev[s].x= v_prev[s].x=dens[s].x= dens_prev[s].x =u[s].x;
       v[s].y = u_prev[s].y = v_prev[s].y= v_prev[s].y=dens[s].y= dens_prev[s].y =u[s].y;
       kL = sizein ;//sobre el borde a la derecha de la esquina
       kU = sizein + 4*N-1;//sobre el borde abajo de la esquina
       u[s].U         = u[s].D         = &u[kU];        
       u[s].L         = u[s].R         = &u[kL];
       v[s].U         = v[s].D         = &v[kU];
       v[s].L         = v[s].R         = &v[kL];
       u_prev[s].U    = u_prev[s].D    = &u_prev[kU];
       u_prev[s].L    = u_prev[s].R    = &u_prev[kL];
       v_prev[s].U    = v_prev[s].D    = &v_prev[kU];
       v_prev[s].L    = v_prev[s].R    = &v_prev[kL];
       dens[s].U      = dens[s].D      = &dens[kU];
       dens[s].L      = dens[s].R      = &dens[kL];
       dens_prev[s].U = dens_prev[s].D = &dens_prev[kU];
       dens_prev[s].L = dens_prev[s].R = &dens_prev[kL];

    
       // esquina arriba derecha
       s = sizein + 4*N + 3;
       u[s].x = N;
       u[s].y = N;
       v[s].x = u_prev[s].x = v_prev[s].x= v_prev[s].x=dens[s].x= dens_prev[s].x =u[s].x;
       v[s].y = u_prev[s].y = v_prev[s].y= v_prev[s].y=dens[s].y= dens_prev[s].y =u[s].y;
       kL =  sizein + N-1;// en el borde a la izquierda de la esquina
       kU =  sizein + 3*N-1;//en el borde abajo de la esquina 
       u[s].U         = u[s].D         = &u[kU];        
       u[s].L         = u[s].R         = &u[kL];
       v[s].U         = v[s].D         = &v[kU];
       v[s].L         = v[s].R         = &v[kL];
       u_prev[s].U    = u_prev[s].D    = &u_prev[kU];
       u_prev[s].L    = u_prev[s].R    = &u_prev[kL];
       v_prev[s].U    = v_prev[s].D    = &v_prev[kU];
       v_prev[s].L    = v_prev[s].R    = &v_prev[kL];
       dens[s].U      = dens[s].D      = &dens[kU];
       dens[s].L      = dens[s].R      = &dens[kL];
       dens_prev[s].U = dens_prev[s].D = &dens_prev[kU];
       dens_prev[s].L = dens_prev[s].R = &dens_prev[kL];
*/

	      
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
        if (max_velocity2 < u[i].val * u[i].val + v[i].val * v[i].val) {
            max_velocity2 = u[i].val * u[i].val + v[i].val * v[i].val;
        }
        if (max_density < d[i].val) {
            max_density = d[i].val;
        }
    }

    for (i = 0; i < size; i++) {
        u[i].val = v[i].val = d[i].val = 0.0f;
    }

    if (max_velocity2 < 0.0000005f) {
        u[sizen/2].val = force * 10.0f;
        v[sizen/2].val = force * 10.0f;
    }
    if (max_density < 1.0f) {
        d[sizen/2].val = source * 10.0f;
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
        for (int i = 0; i <N*N; i++)
        {
                int ii,jj;
                peano_keyto2d(N,i,&ii, &jj); 
                fprintf(fp, "%d %d %e %e %e\n",ii+1,jj+1,dens[i].val,u[i].val,v[i].val);
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

    if (argc != 1 && argc != 7 && argc !=2) {
        fprintf(stderr, "usage : %s N dt diff visc force source\n", argv[0]);
        fprintf(stderr, "where:\n");
        fprintf(stderr, "\t N      : grid resolution\n");
        fprintf(stderr, "\t dt     : time step\n");
        fprintf(stderr, "\t diff   : diffusion rate of the density\n");
        fprintf(stderr, "\t visc   : viscosity of the fluid\n");
        fprintf(stderr, "\t force  : scales the mouse movement that generate a force\n");
        fprintf(stderr, "\t source : amount of density that will be deposited\n");
        exit(1);
    }

    if (argc == 1) {
        N = 128;
        dt = 0.1f;
        diff = 0.0f;
        visc = 0.0f;
        force = 5.0f;
        source = 100.0f;
        fprintf(stderr, "Using defaults : N=%d dt=%g diff=%g visc=%g force = %g source=%g\n",
                N, dt, diff, visc, force, source);
    } 
    if (argc == 2) {
        N = atoi(argv[1]);
        dt = 0.1f;
        diff = 0.0f;
        visc = 0.0f;
        force = 5.0f;
        source = 100.0f;
        fprintf(stderr, "Using defaults : N=%d dt=%g diff=%g visc=%g force = %g source=%g\n",
                N, dt, diff, visc, force, source);
    }else {
        N = atoi(argv[1]);
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

    for (i = 0; i < 100; i++) {
        one_step(i,0);
    }

   free_data();
    exit(0);
}
