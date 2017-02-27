#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include <string>
#include "diffeqsolver.h"
#include "NRutil.h"
#include "process_functions.h"
using namespace std;

void rk4(double *y, double *dydx, int n, double x, double h, double *yout, void (*derivs)(double, double *, double *)) {
    int i;
    double xh,hh,h6;

    double * dym;
    double * dyt;
    double * yt;
    dym = new double [n];
    dyt = new double [n];
    yt = new double [n];

    hh=h*0.5;
    h6=h/6.0;
    xh=x+hh;
    for (i=0;i<n;i++)
        yt[i]=y[i]+hh*dydx[i];
    (*derivs)(xh, yt, dyt);
    for (i=0;i<n;i++)
        yt[i]=y[i]+hh*dyt[i];
    (*derivs)(xh, yt, dym);
    for (i=0;i<n;i++) {
        yt[i]=y[i]+h*dym[i];
        dym[i] += dyt[i];
    }
    (*derivs)(x+h, yt, dyt);
    for (i=0;i<n;i++) yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);

    delete [] dym;
    delete [] dyt;
    delete [] yt;
}

#define PGROW -0.20
#define PSHRNK -0.25
#define FCOR 0.06666666
#define SAFETY 0.9
#define ERRCON 1.89e-4
#define HMAX 10.0

#define MAXSTP 10000 // normally 10000
#define TINY 1.0e-30

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
          ) {
    int i;
    double xsav,hh,h,temp,errmax;

    double * dysav;
    double * ysav;
    double * ytemp;
    dysav = new double [n];
    ysav = new double [n];
    ytemp = new double [n];

    xsav=(*x);
    for (i=0;i<n;i++) {
        ysav[i]=y[i];
        dysav[i]=dydx[i];
    }
    h=htry;
    for (;;) {
        hh=0.5*h;
        rk4(ysav, dysav, n, xsav, hh, ytemp, derivs);
        *x=xsav+hh;
        (*derivs)(*x, ytemp, dydx);
        rk4(ytemp, dydx, n, *x, hh, y, derivs);
        *x=xsav+h;
        if (*x == xsav) {
            cout << "Step size too small in routine RKQC\n";
            break;
        }
        rk4(ysav, dysav, n, xsav, h, ytemp, derivs);

        errmax=0.0;
        for (i=0;i<n;i++) {
            ytemp[i]=y[i]-ytemp[i];
            temp=fabs(ytemp[i]/yscal[i]);
            if (errmax < temp) errmax=temp;
        }
	//cout << "x:" << *x << "\t\tp.a.:" << y[0]*180.0/M_PI << "\t\tV:" << tanh(2*y[1]) << "\t\tbeta_b + delta:" << ( BetaB(*x) + delta(*x) ) * 180.0 / M_PI + 90.0 << "\n";
	//cout << "x:" << *x << "\t\tratio:" << ( Lambda(*x) / Q(*x) ) / ( Lambda(*x) * cos(2.0 * (y[0] - BetaB(*x) - delta(*x))) * sinh(2 * y[1]) ) << "\n";
        errmax /= eps;
        if (errmax <= 1.0)	{
            *hdid=h;
            *hnext=(errmax > ERRCON ?
                    SAFETY*h*exp(PGROW*log(errmax)) : 4.0*h);

            // Making sure that stepsize is not too long

            if (*hnext > HMAX) *hnext = HMAX;
            break;
        }
        h=SAFETY*h*exp(PSHRNK*log(errmax));
    }
    for (i=0;i<n;i++) y[i] += ytemp[i]*FCOR;

    delete [] dysav;
    delete [] ysav;
    delete [] ytemp;
}

#undef PGROW
#undef PSHRNK
#undef FCOR
#undef SAFETY
#undef ERRCON

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
            ) {
    int nstp,i;
    double xsav,x,hnext,hdid,h;

    int kmax = 2, kount = 0;
    double *xp, **yp, dxsav;
    xp = Nvector(0, kmax);
    yp = Nmatrix(0, nvar, 0, kmax);

    double * yscal;
    double * y;
    double * dydx;
    yscal = new double [nvar];
    y = new double [nvar];
    dydx = new double [nvar];

    string path = "dats/test/";
    ofstream plot(path + "betadelta_5.dat");

    x=x1;
    h=(x2 > x1) ? fabs(h1) : -fabs(h1);
    nok = nbad = kount = 0;
    for (i=0;i<nvar;i++) y[i]=ystart[i];
    if (kmax > 0) xsav=x-dxsav*2.0;
    for (nstp=1;nstp<=MAXSTP;nstp++) {
        (*derivs)(x, y, dydx);
        for (i=0;i<nvar;i++)
            yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
        if (kmax > 0) {
            if (fabs(x-xsav) > fabs(dxsav)) {
                if (kount < kmax-1) {
                    xp[++kount]=x;
                    for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
                    xsav=x;
                }
            }
        }
        if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;

        (rkqc)(y, dydx, nvar, &x, h, eps, yscal, &hdid, &hnext, derivs);

        plot << x << " " << BetaB(x) + delta(x) << " " << y[0] << " " << tanh(2 * y[1]) << "\n";

        if (hdid == h) nok++; else nbad++;
        if ((x-x2)*(x2-x1) >= 0.0) {
            for (i=0; i<nvar; i++) {
                ystart[i]=y[i];
            }
            if (kmax) {
                xp[++kount]=x;
                for (i=1; i<=nvar; i++) yp[i][kount]=y[i];
            }
            plot.close();
            return;
        }
        if (fabs(hnext) <= hmin) {
            cout << "Step size too small in ODEINT\n";
            break;
        }
        h=hnext;
    }
    cout << "Too many steps in routine ODEINT\n";
    plot.close();
    delete [] yscal;
    delete [] y;
    delete [] dydx;
}

#undef MAXSTP
#undef TINY
