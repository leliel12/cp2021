//
// solver.h
//

#ifndef SOLVER_H_INCLUDED
#define SOLVER_H_INCLUDED
struct float_link
{
	float val;
	int x,y;
	struct float_link *U;
	struct float_link *D;
	struct float_link *L;
	struct float_link *R;
};

void dens_step(unsigned int n, struct float_link *x, struct float_link *x0, struct float_link *u, struct float_link *v, float diff, float dt);
void vel_step(unsigned int n, struct float_link *u, struct float_link *v, struct float_link *u0, struct float_link *v0, float visc, float dt);
void rot(int n, int *x, int *y, int rx, int ry); 
int peano_2dtokey(int n, int x, int y);
void peano_keyto2d(int n, int d, int *x, int *y);

void init_pointers(struct float_link *ptr, int n);
#endif /* SOLVER_H_INCLUDED */
