#pragma DIFFEQSOLVER

void rk4(double *y, double *dydx, int n, double x, double h, double *yout, void (*derivs)(double, double *, double *));
void rkqc(
          double *y,
          double *dydx,
          int n,
          double *x,
          double htry,
          double eps,
          double *yscal,
          double *hdid,
          double *hnext,
          void (*derivs)(double, double *, double *)
          );

//extern int kmax;
//extern int kount;
//extern double *xp;
//extern double **yp;
//extern double dxsav;

void odeint(
            double ystart[],
            int nvar,
            double x1,
            double x2,
            double eps,
            double h1,
            double hmin,
            int nok,
            int nbad,
            void (*derivs)(double, double *, double *)
            );