#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include "solver_peano.h"

#define SWAP(x0, x)      \
    {                    \
        struct float_link *tmp = x0; \
        x0 = x;          \
        x = tmp;         \
    }

//funciones para ordenar siguiendo la curva de peano-hilbert
//sacada de wikipedia
void rot(int n, int *x, int *y, int rx, int ry) 
{
    if (ry == 0) {
        if (rx == 1) {
            *x = n-1 - *x;
            *y = n-1 - *y;
        }

        //Swap x and y
        int t  = *x;
        *x = *y;
        *y = t;
    }
}

int xy2key(int n, int x, int y)
{
	int sizein = n*n;
	if(x>=0 && x<n) 
	{
		if(y>=0 && y< n)
		{
			return(peano_2dtokey(n,x,y));
		} 
		else if (y==(-1))
		{
	                //borde inferior
			return(sizein + n+ x);
		} 
		else if (y==n)
		{
			//borde superior
			return(sizein + x);
		}
		else
		{
			fprintf(stderr,"Wrong index 1\n");
			exit(-1);
		}

	}
	else 
	{
		if(y>=0 && y < n)
		{
		    if (x==(-1))
		    {
			  //borde izquierdo
	                  return(sizein + 3*n+ y);
		    }
		    else if (x==n) 
		    {
			  //borde derecho
			  return(sizein + 2*n+ y);
		    }
		    else 
		    {
			fprintf(stderr,"Wrong index 2\n");
			exit(-1);
		    }

		} 
		else if (y==(-1))
		{
     		    if (x==(-1))
		    {
			    //esquina abajo izquierda
			    return(sizein + 4*n);
		    } 
		    else if (x==n) 
		    {
			    //esquina abajo derecha
			    return(sizein + 4*n +1);
		    } 
		    else 
		    {
			fprintf(stderr,"Wrong index 3\n");
			exit(-1);
		    }


		} 
		else if (y==n)
		{
		    if (x==(-1))
		    {
			    //esquina arriba izquierda
			    return(sizein + 4*n + 2);
		    } 
		    else if (x==n) 
		    {
			    //esquina arriba derecha
			    return(sizein + 4*n + 3);
		    } 
		    else 
		    {
			fprintf(stderr,"Wrong index 4\n");
			exit(-1);
		    }

		}
		else
		{
			fprintf(stderr,"Wrong index 5 %d %d\n",x,y);
			exit(-1);
		}
	}
	fprintf(stderr,"Wrong xy2key\n");
	exit(-1);
	return(0);
}

int peano_2dtokey(int n, int x, int y) 
{
    int rx, ry, s, d=0;
    for (s=n/2; s>0; s/=2) {
        rx = (x & s) > 0;
        ry = (y & s) > 0;
        d += s * s * ((3 * rx) ^ ry);
        rot(n, &x, &y, rx, ry);
    }
    return d;
}

void peano_keyto2d(int n, int d, int *x, int *y) 
{
    int rx, ry, s, t=d;
    *x = *y = 0;
    for (s=1; s<n; s*=2) {
        rx = 1 & (t/2);
        ry = 1 & (t ^ rx);
        rot(s, x, y, rx, ry);
        *x += s * rx;
        *y += s * ry;
        t /= 4;
    }
}

void init_pointers(struct float_link *ptr,int N)
{

    int k,s;
    int i,j;
    int kU,kD,kL,kR;
    int size   = (N + 2) * (N + 2);
    int sizein = N*N;

    for (i = 0; i <N; i++) 
    {
       for (j = 0; j <N; j++) 
       {
	       k  = xy2key(N,i,j);
	       kU = xy2key(N,i   ,j+1);
	       kD = xy2key(N,i   ,j-1);
	       kR = xy2key(N,i+1 ,j);
	       kL = xy2key(N,i-1 ,j);
               //fprintf(fp, "%d %d %d %d %d %d %d\n",i,j,k,kU,kD,kL,kR);
	       ptr[k].x = i;
	       ptr[k].y = j;
	       ptr[k].U = &ptr[kU];
	       ptr[k].D = &ptr[kD];
	       ptr[k].L = &ptr[kL];
	       ptr[k].R = &ptr[kR];
       }
    }

    //boundaries
    i=0;
    for (int t = 0; t <N; t++) 
    {
	       int i,j;
	       // borde superior 
	       j=N; i=t;
	       k  = xy2key(N,i,j);
	       kD = xy2key(N,i,j-1);
	       kU = kD;//ya estamos arriba apuntar para abajo
	       kR = xy2key(N,i+1,j);
	       kL = xy2key(N,i-1,j);
               //fprintf(fp, "%d %d %d %d %d %d %d\n",i,j,k,kU,kD,kL,kR);
	       ptr[k].x = i;
	       ptr[k].y = j;
	       ptr[k].U = &ptr[kU];
	       ptr[k].D = &ptr[kD];
	       ptr[k].L = &ptr[kL];
	       ptr[k].R = &ptr[kR];

	       // borde inferior j=-1
	       j=-1; i=t;
	       k  = xy2key(N,i,j);
	       kU = xy2key(N,i,j+1);
	       kD = kU;//ya estamos abajo apuntar para arriba
	       kR = xy2key(N,i+1,j);
	       kL = xy2key(N,i-1,j);
               //fprintf(fp, "%d %d %d %d %d %d %d\n",i,j,k,kU,kD,kL,kR);
	       ptr[k].x = i;
	       ptr[k].y = j;
	       ptr[k].U = &ptr[kU];
	       ptr[k].D = &ptr[kD];
	       ptr[k].L = &ptr[kL];
	       ptr[k].R = &ptr[kR];

   	       // borde derecho i=N
	       i=N; j=t;
	       k  = xy2key(N,i,j);
	       kU = xy2key(N,i,j+1);
	       kD = xy2key(N,i,j-1);
	       kL = xy2key(N,i-1,j);
	       kR = kL;//ya estamos a la derecha apuntar a la izquierda
               //fprintf(fp, "%d %d %d %d %d %d %d\n",i,j,k,kU,kD,kL,kR);
	       ptr[k].x = i;
	       ptr[k].y = j;
	       ptr[k].U = &ptr[kU];
	       ptr[k].D = &ptr[kD];
	       ptr[k].L = &ptr[kL];
	       ptr[k].R = &ptr[kR];


   	       // borde izquierdo i=-1
	       i=-1; j =t;
	       k  = xy2key(N,i,j);
	       kU = xy2key(N,i,j+1);
	       kD = xy2key(N,i,j-1);
	       kR = xy2key(N,i+1,j);
	       kL = kR;//ya estamos a la izquierda apuntar a la derecha
               //fprintf(fp, "%d %d %d %d %d %d %d\n",i,j,k,kU,kD,kL,kR);
	       ptr[k].x = i;
	       ptr[k].y = j;
	       ptr[k].U = &ptr[kU];
	       ptr[k].D = &ptr[kD];
	       ptr[k].L = &ptr[kL];
	       ptr[k].R = &ptr[kR];

    }
    // esquinas

    // esquina abajo izquierda
    i=-1;j=-1;
    k  = xy2key(N,i,j);
    kU = xy2key(N,i,j+1);
    kD = kU;
    kR = xy2key(N,i+1,j);
    kL = kR;
    //fprintf(fp, "%d %d %d %d %d %d %d\n",i,j,k,kU,kD,kL,kR);
    ptr[k].x = i;
    ptr[k].y = j;
    ptr[k].U = &ptr[kU];
    ptr[k].D = &ptr[kD];
    ptr[k].L = &ptr[kL];
    ptr[k].R = &ptr[kR];
   
    // esquina abajo derecha
    i=N;j=-1;
    k  = xy2key(N,i,j);
    kU = xy2key(N,i,j+1);
    kD = kU;
    kL = xy2key(N,i-1,j);
    kR = kL;
    //fprintf(fp, "%d %d %d %d %d %d %d\n",i,j,k,kU,kD,kL,kR);
    ptr[k].x = i;
    ptr[k].y = j;
    ptr[k].U = &ptr[kU];
    ptr[k].D = &ptr[kD];
    ptr[k].L = &ptr[kL];
    ptr[k].R = &ptr[kR];
   
       
    // esquina arriba izquierda
    i=-1;j=N;
    k  = xy2key(N,i,j);
    kD = xy2key(N,i,j-1);
    kU = kD;
    kR = xy2key(N,i+1,j);
    kL = kR;
    //fprintf(fp, "%d %d %d %d %d %d %d\n",i,j,k,kU,kD,kL,kR);
    ptr[k].x = i;
    ptr[k].y = j;
    ptr[k].U = &ptr[kU];
    ptr[k].D = &ptr[kD];
    ptr[k].L = &ptr[kL];
    ptr[k].R = &ptr[kR];
       
    // esquina arriba derecha
    i=N; j=N;
	       
    k  = xy2key(N,i,j);
    kD = xy2key(N,i,j-1);
    kU = kD;
    kL = xy2key(N, i-1, j   );
    kR = kL;
    //fprintf(fp, "%d %d %d %d %d %d %d\n",i,j,k,kU,kD,kL,kR);
    ptr[k].x = i;
    ptr[k].y = j;
    ptr[k].U = &ptr[kU];
    ptr[k].D = &ptr[kD];
    ptr[k].L = &ptr[kL];
    ptr[k].R = &ptr[kR];


}

typedef enum 
{ NONE = 0,
  VERTICAL = 1,
  HORIZONTAL = 2 }
  boundary;

static void add_source(unsigned int n,struct float_link *x, const struct float_link *s, float dt)
{
    unsigned int size = (n + 2) * (n + 2);
    for (unsigned int i = 0; i < size; i++) {
        x[i].val += dt * s[i].val;
    }
}

static void set_bnd(unsigned int n, boundary b, struct float_link *x)
{
    int size   = (n + 2) * (n + 2);
    int sizein = n*n;
    int f1,f2;

    if(b ==VERTICAL)
    {
	    f1 = 1;
	    f2 =-1;
    } 
    else if (b==HORIZONTAL)
    {
	    f1 = -1;
	    f2 =  1;
    }
    else
    {
	    f1 = 1;
	    f2 = 1;
    }
    // borde superior e inferior (HORIZONTAL)
    for (int t = sizein; t < 2*n+sizein; t++) 
	    x[t].val = f1*((x[t].U)-> val);
	    
    // borde derecho e izquierdo (VERTICAL)
    for (int t = 2*n+sizein; t <4*n+sizein; t++) 
	    x[t].val = f2*(((x[t].R)-> val));

    // esquinas
    for(int t=size-4; t<size; t++)
	    x[t].val = 0.5f*(((x[t].D)-> val)+((x[t].R)-> val));
	
}

static void lin_solve(unsigned int n, boundary b, struct float_link *x, const struct float_link* x0, float a, float c)
{
    int size = (n + 2) * (n + 2);
    int sizein = n*n;
    struct float_link * xk;

    //que siempre sea par para que el swap xk no quede invertido
    for (int k = 0; k < 20; k++) 
    {
        for (int i = 0; i < sizein; i=i+2) 
	{
	      x[i].val = (x0[i].val + a * ((x[i].L)->val + (x[i].R)->val + (x[i].D)->val + (x[i].U)->val)) / c;
        }

        for (int i = 1; i < sizein; i=i+2) 
	{
	      x[i].val = (x0[i].val + a * ((x[i].L)->val + (x[i].R)->val + (x[i].D)->val + (x[i].U)->val)) / c;
        }
        set_bnd(n, b, x);
    }
}

static void diffuse(unsigned int n, boundary b, struct float_link *x, const struct float_link *x0, float diff, float dt)
{
    float a = dt * diff * n * n;
    lin_solve(n, b, x, x0, a, 1 + 4 * a);
}

static void advect(unsigned int n, boundary b, struct float_link *d, struct float_link * d0, const struct float_link* u, const struct float_link* v, float dt)
{
    int i0, j0;
    int i1, j1;
    unsigned int size = (n + 2) * (n + 2);
    unsigned int sizein =  n*n;
    int sx;
    int sy;
    struct float_link *dd0;
    float x, y, s0, t0, s1, t1;

    float dt0 = dt * n;
    for (unsigned int i = 0; i < sizein; i++) 
    {
	    x = (d0[i].x+1) - dt0 * u[i].val;
            y = (d0[i].y+1) - dt0 * v[i].val;
            if (x < 0.5f) {
                x = 0.5f;
            } else if (x > n + 0.5f) {
                x = n + 0.5f;
            }
            if (y < 0.5f) {
                y = 0.5f;
            } else if (y > n + 0.5f) {
                y = n + 0.5f;
            }
            i0 = (int)x;
            i1 = i0 + 1;
	    j0 = (int)y;
            j1 = j0 + 1;
            s1 = x - i0;
            s0 = 1 - s1;
            t1 = y - j0;
            t0 = 1 - t1;
	    dd0=&d0[xy2key(n,i0-1,j0-1)];

	    d[i].val = s0 * (t0 * dd0->val     + t1 * (dd0->U)->val) 
		     + s1 * (t0 *(dd0->R)->val + t1 * ((dd0->R)->U)->val);
    }

    set_bnd(n, b, d);
}

static void project(unsigned int n,struct float_link * u, struct float_link *v, struct float_link *p, struct float_link * div)
{
    unsigned int sizein =  n*n;
    for (unsigned int i = 0; i < sizein; i++) 
    {
	    div[i].val = -0.5f * ((u[i].R)->val - (u[i].L)->val + (v[i].U)->val - (v[i].D)->val) / n;
            p[i].val = 0.0f;
    }
    set_bnd(n, NONE, div);
    set_bnd(n, NONE, p);

    lin_solve(n, NONE, p, div, 1, 4);

    for (unsigned int i = 0; i < sizein; i++) 
    {
	    u[i].val -= 0.5f * n * ((p[i].R)->val - (p[i].L)->val);
	    v[i].val -= 0.5f * n * ((p[i].U)->val - (p[i].D)->val);
    }
    set_bnd(n, VERTICAL, u);
    set_bnd(n, HORIZONTAL, v);

}

void dens_step(unsigned int n, struct float_link *x, struct float_link *x0, struct float_link *u, struct float_link *v, float diff, float dt)
{
    add_source(n, x, x0, dt);
    SWAP(x0, x);
    diffuse(n, NONE, x, x0, diff, dt);
    SWAP(x0, x);
    advect(n, NONE, x, x0, u, v, dt);

}

void vel_step(unsigned int n, struct float_link *u, struct float_link *v, struct float_link *u0, struct float_link *v0, float visc, float dt)
{
    add_source(n, u, u0, dt);
    add_source(n, v, v0, dt);
    SWAP(u0, u);
    diffuse(n, VERTICAL, u, u0, visc, dt);
    SWAP(v0, v);
    diffuse(n, HORIZONTAL, v, v0, visc, dt);
    project(n, u, v, u0, v0);
    SWAP(u0, u);
    SWAP(v0, v);
    advect(n, VERTICAL, u, u0, u0, v0, dt);
    advect(n, HORIZONTAL, v, v0, u0, v0, dt);
    project(n, u, v, u0, v0);
}
