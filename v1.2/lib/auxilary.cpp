#include <iostream>
#include <math.h>
#include <vector>
#include "pulsar_parameters.h"
#include "constants.h"
#include "functions.h"
#include "auxilary.h"
using namespace std;

double theta_em;
double phi_em;

// finding initial point of generation (along field line)

vector <double> b0 (double th, double ph, double PHI0) {
    vector <double> n0(3);
    n0[0] = sin(th) * cos(ph);
    n0[1] = sin(th) * sin(ph);
    n0[2] = cos(th);

    vector <double> m(3);
    m[0] = sin(alpha) * cos(PHI0);
    m[1] = sin(alpha) * sin(PHI0);
    m[2] = cos(alpha);

    return NORMALIZE(SUM(TIMES(3.0 * SCALAR(m, n0), n0), TIMES(-1.0, m)));
}
double func1 (double x, double y, double PHI0) {
    vector <double> o(3);
    o[0] = sin(dzeta);
    o[1] = 0.0;
    o[2] = cos(dzeta);
    return CROSS(b0(x, y, PHI0), o)[0];
}
double func2 (double x, double y, double PHI0) {
    vector <double> o(3);
    o[0] = sin(dzeta);
    o[1] = 0.0;
    o[2] = cos(dzeta);
    return CROSS(b0(x, y, PHI0), o)[1];
}
double DX (double (*func)(double, double, double), double x, double y, double PHI0) {
    double h = 0.00001;
    double fm2 = func(x - 2 * h, y, PHI0);
    double fp2 = func(x + 2 * h, y, PHI0);
    double fm1 = func(x - h, y, PHI0);
    double fp1 = func(x + h, y, PHI0);
    return (fm2 - 8 * fm1 + 8 * fp1 - fp2) / (12 * h);
}
double DY (double (*func)(double, double, double), double x, double y, double PHI0) {
    double h = 0.00001;
    double fm2 = func(x, y - 2 * h, PHI0);
    double fp2 = func(x, y + 2 * h, PHI0);
    double fm1 = func(x, y - h, PHI0);
    double fp1 = func(x, y + h, PHI0);
    return (fm2 - 8 * fm1 + 8 * fp1 - fp2) / (12 * h);
}

void findInitPoints (double PHI0) {
    double X = alpha, Y = PHI0;
    for(int i = 0; i < 100; i ++) {
        double f1x = DX(func1, X, Y, PHI0);
        double f2x = DX(func2, X, Y, PHI0);
        double f1y = DY(func1, X, Y, PHI0);
        double f2y = DY(func2, X, Y, PHI0);
        double f1 = func1(X, Y, PHI0);
        double f2 = func2(X, Y, PHI0);
        double dX = (f1y * f2 - f1 * f2y) / (f1x * f2y - f1y * f2x);
        double dY = (f1x * f2 - f1 * f2x) / (f1y * f2x - f2y * f1x);
        X += dX;
        Y += dY;
    }
    theta_em = X;
    phi_em = Y;
}
