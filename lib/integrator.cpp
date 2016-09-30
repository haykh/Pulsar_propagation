#include "integrator.h"
//#include <math.h>

double integrate (double (*func)(double), double a, double b) {
	double r = a, dr = (b - a) / 500.0;
	double s = 0.0;
	while (r < b - dr) {
		s += 0.5 * dr * (func(r) + func(r + dr));
		r += dr;
	}
	return s;
}

//double trapzd(double (*func)(double, double), double a, double b, int n, double PHI0) {
//    double x,tnm,sum,del;
//    static double s;
//    int it,j;
//    
//    if (n == 1) {
//        return (s=0.5*(b-a)*(func(a, PHI0)+func(b, PHI0)));
//    } else {
//        for (it=1,j=1;j<n-1;j++) it <<= 1;
//        tnm=it;
//        del=(b-a)/tnm;
//        x=a+0.5*del;
//        for (sum=0.0,j=1;j<=it;j++,x+=del) sum += func(x, PHI0);
//        s=0.5*(s+(b-a)*sum/tnm);
//        return s;
//    }
//}
//
//#define EPS 1.0e-6
//#define JMAX 200
//
//double integrate(double (*func)(double, double), double a, double b, double PHI0) {
//    double trapzd(double (*func)(double, double), double a, double b, int n, double PHI0);
//    int j;
//    double s,st,ost,os;
//
//    ost = os = -1.0e30;
//    for (j = 1; j <= JMAX; j ++) {
//        st=trapzd(func, a, b, j, PHI0);
//        s=(4.0*st-ost)/3.0;
//        if (fabs(s-os) < EPS*fabs(os)) return s;
//        os=s;
//        ost=st;
//    }
//    return 0.0;
//}
//#undef EPS
//#undef JMAX