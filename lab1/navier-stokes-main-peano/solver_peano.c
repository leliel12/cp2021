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

    for (i = 1; i <=N; i++) 
    {
       for (j = 1; j <=N; j++) 
       {
	       k  = peano_2dtokey(N+2,i,j);
	       kU = peano_2dtokey(N+2,i   ,j+1);
	       kD = peano_2dtokey(N+2,i   ,j-1);
	       kR = peano_2dtokey(N+2,i+1 ,j);
	       kL = peano_2dtokey(N+2,i-1 ,j);
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
    for (int t = 1; t <=N; t++) 
    {
	       int i,j;
	       // borde superior 
	       j=N+1; i=t;
	       k  = peano_2dtokey(N+2,i,j);
	       kD = peano_2dtokey(N+2,i,j-1);
	       kU = kD;//ya estamos arriba apuntar para abajo
	       kR = peano_2dtokey(N+2,i+1,j);
	       kL = peano_2dtokey(N+2,i-1,j);
               //fprintf(fp, "%d %d %d %d %d %d %d\n",i,j,k,kU,kD,kL,kR);
	       ptr[k].x = i;
	       ptr[k].y = j;
	       ptr[k].U = &ptr[kU];
	       ptr[k].D = &ptr[kD];
	       ptr[k].L = &ptr[kL];
	       ptr[k].R = &ptr[kR];

	       // borde inferior j=0
	       j=0; i=t;
	       k  = peano_2dtokey(N+2,i,j);
	       kU = peano_2dtokey(N+2,i,j+1);
	       kD = kU;//ya estamos abajo apuntar para arriba
	       kR = peano_2dtokey(N+2,i+1,j);
	       kL = peano_2dtokey(N+2,i-1,j);
               //fprintf(fp, "%d %d %d %d %d %d %d\n",i,j,k,kU,kD,kL,kR);
	       ptr[k].x = i;
	       ptr[k].y = j;
	       ptr[k].U = &ptr[kU];
	       ptr[k].D = &ptr[kD];
	       ptr[k].L = &ptr[kL];
	       ptr[k].R = &ptr[kR];

   	       // borde derecho i=N+1
	       i=N+1; j=t;
	       k  = peano_2dtokey(N+2,i,j);
	       kU = peano_2dtokey(N+2,i,j+1);
	       kD = peano_2dtokey(N+2,i,j-1);
	       kL = peano_2dtokey(N+2,i-1,j);
	       kR = kL;//ya estamos a la derecha apuntar a la izquierda
               //fprintf(fp, "%d %d %d %d %d %d %d\n",i,j,k,kU,kD,kL,kR);
	       ptr[k].x = i;
	       ptr[k].y = j;
	       ptr[k].U = &ptr[kU];
	       ptr[k].D = &ptr[kD];
	       ptr[k].L = &ptr[kL];
	       ptr[k].R = &ptr[kR];


   	       // borde izquierdo i=0
	       i=0; j =t;
	       k  = peano_2dtokey(N+2,i,j);
	       kU = peano_2dtokey(N+2,i,j+1);
	       kD = peano_2dtokey(N+2,i,j-1);
	       kR = peano_2dtokey(N+2,i+1,j);
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
    i=0;j=0;
    k  = peano_2dtokey(N+2,i,j);
    kU = peano_2dtokey(N+2,i,j+1);
    kD = kU;
    kR = peano_2dtokey(N+2,i+1,j);
    kL = kR;
    //fprintf(fp, "%d %d %d %d %d %d %d\n",i,j,k,kU,kD,kL,kR);
    ptr[k].x = i;
    ptr[k].y = j;
    ptr[k].U = &ptr[kU];
    ptr[k].D = &ptr[kD];
    ptr[k].L = &ptr[kL];
    ptr[k].R = &ptr[kR];
   
    // esquina abajo derecha
    i=N+1;j=0;
    k  = peano_2dtokey(N+2,i,j);
    kU = peano_2dtokey(N+2,i,j+1);
    kD = kU;
    kL = peano_2dtokey(N+2,i-1,j);
    kR = kL;
    //fprintf(fp, "%d %d %d %d %d %d %d\n",i,j,k,kU,kD,kL,kR);
    ptr[k].x = i;
    ptr[k].y = j;
    ptr[k].U = &ptr[kU];
    ptr[k].D = &ptr[kD];
    ptr[k].L = &ptr[kL];
    ptr[k].R = &ptr[kR];
   
       
    // esquina arriba izquierda
    i=0;j=N+1;
    k  = peano_2dtokey(N+2,i,j);
    kD = peano_2dtokey(N+2,i,j-1);
    kU = kD;
    kR = peano_2dtokey(N+2,i+1,j);
    kL = kR;
    //fprintf(fp, "%d %d %d %d %d %d %d\n",i,j,k,kU,kD,kL,kR);
    ptr[k].x = i;
    ptr[k].y = j;
    ptr[k].U = &ptr[kU];
    ptr[k].D = &ptr[kD];
    ptr[k].L = &ptr[kL];
    ptr[k].R = &ptr[kR];
       
    // esquina arriba derecha
    i=N+1; j=N+1;
	       
    k  = peano_2dtokey(N+2,i,j);
    kD = peano_2dtokey(N+2,i,j-1);
    kU = kD;
    kL = peano_2dtokey(N+2, i-1, j   );
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
        x[i].val[0] += dt * s[i].val[0];
    }
}

static void lin_solve(unsigned int n, boundary b, struct float_link *x, const struct float_link* x0, float a, float c)
{
    int size = (n + 2) * (n + 2);
    struct float_link * xk;
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

    for (int k = 0; k < 21; k++) 
    {
	int prev= 0 ;//k%2;
	int step= 0 ;//!prev;
	//indices pares
        for (int s = 0; s < size; s=s+2) 
	{
		int i =x[s].x;
		int j =x[s].y;
		if(i>0 && i<=n)
		{
		   if(j>0 && j<=n)
		   {
	                x[s].val[step] = (x0[s].val[0] + a * ((x[s].L)->val[prev] + (x[s].R)->val[prev] + (x[s].D)->val[prev] + (x[s].U)->val[prev])) / c;
		   }
		} 
        }

        for (int s = 1; s < size; s=s+2) 
	{
		int i =x[s].x;
		int j =x[s].y;
		if(i>0 && i<=n)
		{
		   if(j>0 && j<=n)
		   {
	                x[s].val[step] = (x0[s].val[0] + a * ((x[s].L)->val[prev] + (x[s].R)->val[prev] + (x[s].D)->val[prev] + (x[s].U)->val[prev])) / c;
		   }
		} 
        }
	//set bounds
        for (int s = 0; s < size; s++) 
	{
		int i =x[s].x;
		int j =x[s].y;
		if(i>0 && i<=n)
		{
		   if(!(j>0 && j<=n))
		   {
			x[s].val[step] = f1*((x[s].U)-> val[step]);
		   }
		} 
		else
		{
		   if(j>0 && j<=n)//bordes verticales
		   {
	                x[s].val[step] = f2*(((x[s].R)-> val[step]));
		   } 
		   else  //esquinas
		   { 
			x[s].val[step] = 0.5f*(((x[s].D)-> val[step])+((x[s].R)-> val[step]));
		   }
		}
        }
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
    int sx;
    int sy;
    struct float_link *dd0;
    float x, y, s0, t0, s1, t1;

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


    float dt0 = dt * n;
    for (unsigned int s = 0; s < size; s++) 
    {
	    int i = u[s].x;
	    int j = u[s].y;

	    if(i>0 && i<=n)
	    {
	       if(j>0 && j<=n)
	       {
	            x = i - dt0 * u[s].val[0];
                    y = j - dt0 * v[s].val[0];
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
	            dd0=&d0[peano_2dtokey(n+2,i0,j0)];

	            d[s].val[0] = s0 * (t0 * dd0->val[0]     + t1 * (dd0->U)->val[0]) 
	                     + s1 * (t0 *(dd0->R)->val[0] + t1 * ((dd0->R)->U)->val[0]);

	       }
	    } 
    }

    //set boundary
    for (unsigned int s = 0; s < size; s++) 
    {
	    int i = u[s].x;
	    int j = u[s].y;

	    if(i>0 && i<=n)
	    {
	       if(!(j>0 && j<=n)) //bordes horizontales
	       {
	    	d[s].val[0] = f1*((d[s].U)-> val[0]);
	       }
	    } 
	    else
	    {
	       if(j>0 && j<=n)//bordes verticales
	       {
	            d[s].val[0] = f2*(((d[s].R)-> val[0]));
	       } 
	       else  //esquinas
	       { 
	    	d[s].val[0] = 0.5f*(((d[s].D)-> val[0])+((d[s].R)-> val[0]));
	       }
	    }

    }

}

static void project(unsigned int n,struct float_link * u, struct float_link *v, struct float_link *p, struct float_link * div)
{
    int size = (n+2)*(n+2);

    for (unsigned int s = 0; s < size; s++) 
    {
	    int i = div[s].x;
	    int j = div[s].y;

	    if(i>0 && i<=n)
	    {
	       if(j>0 && j<=n)
	       {
	           div[s].val[0] = -0.5f * ((u[s].R)->val[0] - (u[s].L)->val[0] + (v[s].U)->val[0] - (v[s].D)->val[0]) / n;
                   p[s].val[0] = 0.0f;

	       }
	    } 
    }

    //set boundaries
    for (unsigned int s = 0; s < size; s++) 
    {
	    int i = div[s].x;
	    int j = div[s].y;

	    if(i>0 && i<=n)
	    {
	       if(!(j>0 && j<=n)) //bordes horizontales
	       {
	    	div[s].val[0] = (div[s].U)-> val[0];
	    	p[s].val[0] = (p[s].U)-> val[0];
	       }
	    } 
	    else
	    {
	       if(j>0 && j<=n)//bordes verticales
	       {
	            div[s].val[0] = (div[s].R)-> val[0];
	            p[s].val[0] = (p[s].R)-> val[0];
	       } 
	       else  //esquinas
	       { 
	    	    div[s].val[0] = 0.5f*(((div[s].D)-> val[0])+((div[s].R)-> val[0]));
	    	    p[s].val[0]   = 0.5f*(((p[s].D)-> val[0])+((p[s].R)-> val[0]));
	       }
	    }

    }

    lin_solve(n, NONE, p, div, 1, 4);

    for (unsigned int s = 0; s < size; s++) 
    {
	    int i = u[s].x;
	    int j = u[s].y;

	    if(i>0 && i<=n)
	    {
	       if(j>0 && j<=n)
	       {
	           u[s].val[0] -= 0.5f * n * ((p[s].R)->val[0] - (p[s].L)->val[0]);
	           v[s].val[0] -= 0.5f * n * ((p[s].U)->val[0] - (p[s].D)->val[0]);
	       }
	    } 
    }

    //set boundaries
    for (unsigned int s = 0; s < size; s++) 
    {
	    int i = u[s].x;
	    int j = u[s].y;

	    if(i>0 && i<=n)
	    {
	       if(!(j>0 && j<=n))
	       {
	    	  u[s].val[0] = (u[s].U)-> val[0];
	    	  v[s].val[0] = (-((v[s].U)-> val[0]));
	       }
	    } 
	    else
	    {
	       if(j>0 && j<=n)//bordes verticales
	       {
	            u[s].val[0] = (-( (u[s].R)-> val[0]));
	            v[s].val[0] = (v[s].R)-> val[0];
	       } 
	       else  //esquinas
	       { 
	    	    u[s].val[0] = 0.5f*(((u[s].D)-> val[0])+((u[s].R)-> val[0]));
	    	    v[s].val[0]   = 0.5f*(((v[s].D)-> val[0])+((v[s].R)-> val[0]));
	       }
	    }
    }
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
