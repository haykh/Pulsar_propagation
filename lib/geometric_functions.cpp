#include <math.h>
#include <iostream>
#include <vector>
#include "pulsar_parameters.h"
#include "constants.h"
#include "functions.h"
#include "geometric_functions.h"
#include "physical_functions.h"
using namespace std;

// most of functions should be local

double sgn (double value) {
    if (value >= 0.0) {
        return 1.0;
    } else {
        return -1.0;
    }
}

vector <double> vMoment (double R) {
    vector <double> mvec(3);
    mvec[0] = sin(alpha) * cos(PHI0 + R / RLC);
    mvec[1] = sin(alpha) * sin(PHI0 + R / RLC);
    mvec[2] = cos(alpha);
    return mvec;
} // Moment vector
vector <double> vR (double R) {
    vector <double> n0(3);
    n0[0] = sin(theta_em) * cos(phi_em);
    n0[1] = sin(theta_em) * sin(phi_em);
    n0[2] = cos(theta_em);

    vector <double> o(3);
    o[0] = sin(dzeta);
    o[1] = 0.0;
    o[2] = cos(dzeta);
    return SUM(TIMES(R_em, n0), TIMES(R, o));
} // Propagation radius vector
double psi_m (double R) {
    return ANGLE(vR(R), vMoment(R));
}

vector <double> vBdipole (double R) {
    vector <double> m;
    vector <double> n;
    m = vMoment (R);
    n = NORMALIZE(vR(R));
    return SUM(TIMES(3.0 * SCALAR(m, n), n), TIMES(-1.0, m));
}
vector <double> vBsplit (double R) {
    double Rr = NORM(vR(R));

    double costh = SCALAR(NORMALIZE(vR(R)), NORMALIZE(vOmega));
    double sinth = sqrt(1 - costh * costh);
    double cosphi = vR(R)[0] / Rr;
    double sinphi = vR(R)[1] / Rr;
    double phi = acos (cosphi);

    double psi1 = costh * cos(alpha) + sinth * sin(alpha) * cos(phi - PHI0 + Rr / RLC);

    double Br = (fr / (Rr * Rr * RLC)) * tanh(psi1 / 0.1);
    double Bphi = -(fphi * fr * sinth / (Rr * RLC * RLC)) * tanh(psi1 / 0.1);

    Br *= pow(Rr, 3);
    Bphi *= pow(Rr, 3);

    vector <double> temp(3);
    temp[0] = Br * sinth * cosphi - Bphi * sinphi;
    temp[1] = Br * sinth * sinphi + Bphi * cosphi;
    temp[2] = Br * costh;
    return temp;
}
vector <double> vBtor (double R) {
    double Rr = NORM(vR(R));

    double rho = Rr * sin(psi_m(R));
    double R0 = Rr * sqrt(Rr / RLC);

    vector <double> temp(3);
    if (rho > R0) {
        temp[0] = 0.0;
        temp[1] = 0.0;
        temp[2] = 0.0;
        return temp;
    }

    vector <double> n;
    vector <double> m;
    vector <double> O;
    n = NORMALIZE(vR(R));
    m = NORMALIZE(vMoment(R));
    O = NORMALIZE(vOmega);

    vector <double> nPerp;
    vector <double> OPerp;
    nPerp = NORMALIZE(SUM(n, TIMES(-SCALAR(n, m), m)));
    OPerp = NORMALIZE(SUM(O, TIMES(-SCALAR(O, m), m)));

    double varphi = PI + sgn(SCALAR(CROSS(OPerp, nPerp), m)) * ANGLE(OPerp, nPerp);

    double COEFF = (-3.0 / 4.0) * sin(alpha);
    double Brho = (COEFF * sin(varphi) / Rr) * (rho * rho - R0 * R0);
    double Bvarphi = (COEFF * cos(varphi) / Rr) * (3.0 * rho * rho - R0 * R0);

    double Bx1 = Brho * sin(varphi) + Bvarphi * cos(varphi);
    double By1 = - Brho * cos(varphi) + Bvarphi * sin(varphi);

    temp[0] = -Bx1 * sin(PHI0 + R / RLC) - By1 * cos(alpha) * cos(PHI0 + R / RLC);
    temp[1] = Bx1 * cos(PHI0 + R / RLC) - By1 * cos(alpha) * sin(PHI0 + R / RLC);
    temp[2] = By1 * sin(alpha);

    return TIMES(BMULT, temp);
}

vector <double> vB (double R) {
    return SUM(SUM(vBsplit(R), vBdipole(R)), vBtor(R));
}
vector <double> vb (double R) {
    return NORMALIZE(vB(R));
}

double theta_kb (double R) {
    vector <double> o(3);
    o[0] = sin(dzeta);
    o[1] = 0.0;
    o[2] = cos(dzeta);
    return ANGLE(vB(R), o);
}

vector <double> vBetaR (double R) {
    return TIMES(R_star / c, CROSS(vOmega, vR(R)));
}
vector <double> vUdr (double R) {
    vector <double> o(3);
    o[0] = sin(dzeta);
    o[1] = 0.0;
    o[2] = cos(dzeta);

    vector <double> vn;
    vector <double> vm;
    vn = NORMALIZE(SUM(o, TIMES(-SCALAR(o, vb(R)), vb(R))));
    vm = NORMALIZE(CROSS(vb(R), vn));

    vector <double> temp(3);
    temp[0] = SCALAR(vBetaR(R), vn);
    temp[1] = SCALAR(vBetaR(R), vm);
    temp[2] = sqrt(1 - pow(temp[0], 2) - pow(temp[1], 2));
    return temp;
}
double delta (double R) {
    double vx = vUdr(R) [0];
    double vy = vUdr(R) [1];
    double sinth = sin(theta_kb(R));
    double costh = cos(theta_kb(R));
    double sign = sgn (- vy * costh / sqrt(pow(sinth - vx, 2) + pow(costh * vy , 2)));

//    double addit = 0.0;
//    if (sgn(vUdr(10.0)[1] * vy) == -1.0) { // don't even ask what this is
//        addit = -2.0 * PI;
//    }
    return sign * acos((sinth - vx) / sqrt(pow(sinth - vx, 2) + pow(costh * vy , 2)));
}
double BetaB (double R) {
    vector <double> o(3);
    o[0] = sin(dzeta);
    o[1] = 0.0;
    o[2] = cos(dzeta);

    vector <double> XX;
    vector <double> YY;
    XX = NORMALIZE(SUM(vOmega, TIMES(-SCALAR(o, vOmega), o)));
    YY = CROSS(o, XX);

    double bx = SCALAR(XX, vB(R));
    double by = SCALAR(YY, vB(R));
    double bb = sqrt (bx * bx + by * by);

//    return acos (bx / bb) * sgn (by);
    return atan(by / bx);
}

double gFunc (double f, double R) {
    double theta = ANGLE(vR(R), vOmega);
    double dtheta = 5.0 * PI / 180.0;
    double gap;
    if (alpha_deg > 80)
      gap = (1 - exp(-pow(PI / 2.0 - theta, 2) / (2.0 * dtheta * dtheta)));
    else
      gap = 1.0;
    // double rperp = sin(psi_m(R)) * Rr / (Rr * sqrt(Rr / RLC));
    return (pow(f, 2.5) * exp(-f * f) / (pow(f, 2.5) + pow(f0, 2.5))) * gap;
}
double Ne (double R) {
    double f = pow(sin(psi_m(R)), 2) * RLC / NORM(vR(R));
    double nGJ = SCALAR(vOmega, vB(R)) * (B0 / pow(NORM(vR(R)), 3)) / (2 * PI * c * e);
    return lambda * gFunc (f, R) * nGJ;
}

double omegaB (double R) {
    return -e * NORM(vB(R)) * (B0 / pow(NORM(vR(R)), 3)) / (me * c);
}
double omegaW (double R) {
    double vx = vUdr(R) [0];
    double vz = vUdr(R) [2];
    double sinth = sin(theta_kb(R));
    double costh = cos(theta_kb(R));
    return omega * (1 - sinth * vx - costh * vz);
}
double omegaP (double R) {
    return sqrt(4 * PI * e * e * fabs(Ne(R)) / me);
}

double Q (double R) {
    double vx = vUdr(R) [0];
    double vy = vUdr(R) [1];
    double vz = vUdr(R) [2];
    double sinth = sin(theta_kb(R));
    double costh = cos(theta_kb(R));

    return lambda * omegaB(R) * omega * (pow(sinth - vx, 2) + pow(vy * costh, 2)) / (2 * pow(gamma0, 3) * pow(omegaW(R), 2) * (costh * (1 - vx * vx - vy * vy) - vz * (1.0 - sinth * vx)));
}

double fDist (double gamma) {
    return ((6 * gamma0) / (pow(2.0, 1.0/6.0) * PI)) * (pow(gamma, 4) / (2 * pow(gamma, 6) + pow(gamma0, 6)));
}

double gammaU (double R) {
    double vx = vUdr(R) [0];
    double vy = vUdr(R) [1];
    return pow(1 - vx * vx - vy * vy, -0.5);
}
double A (double R) {
    return pow(gammaU(R) * omegaW(R) / omegaB(R), 2);
}

double INTEGRAL (double gamma, double R) {
    double cA = A(R);
    return -(pow(2, 2.0 / 3.0)*(-2*sqrt(3)*atan((2*pow(2, 1.0 / 3.0)*pow(gamma,2) - pow(gamma0,2))/
            (sqrt(3)*pow(gamma0,2))) - 2*log(pow(2, 1.0 / 3.0)*pow(gamma,2) + pow(gamma0,2)) +
         log(pow(2, 2.0 / 3.0)*pow(gamma,4) - pow(2, 1.0 / 3.0)*pow(gamma,2)*pow(gamma0,2) +
           pow(gamma0,4))) - pow(2, 1.0 / 3.0)*pow(gamma0,2)*cA*
       (2*sqrt(3)*atan((2*pow(2, 1.0 / 3.0)*pow(gamma,2) - pow(gamma0,2))/(sqrt(3)*pow(gamma0,2))) -
         2*log(pow(2, 1.0 / 3.0)*pow(gamma,2) + pow(gamma0,2)) +
         log(pow(2, 2.0 / 3.0)*pow(gamma,4) - pow(2, 1.0 / 3.0)*pow(gamma,2)*pow(gamma0,2) +
           pow(gamma0,4))) - 2*pow(gamma0,4)*pow(cA,2)*
       (log(2*pow(gamma,6) + pow(gamma0,6)) - 3*log(fabs((1 - gamma*sqrt(cA))*(1 + gamma*sqrt(cA))))))/
   (2.*pow(2, 1.0 / 6.0)*PI*pow(gamma0,3)*(2 + pow(gamma0,6)*pow(cA,3)));
}
double Lambda (double R) {
    double vx = vUdr(R) [0];
    double vy = vUdr(R) [1];
    double sinth = sin(theta_kb(R));
    double costh = cos(theta_kb(R));
    double avrg = INTEGRAL(1000000.0, R) - INTEGRAL(0.0, R);
    return (-1.0 / 2.0) * pow(omegaP(R) * gammaU(R) / omegaW(R), 2) * avrg * (pow(sinth - vx, 2) + pow(vy * costh, 2));
}

double dtau (double R) {
    return pow(omegaP(R), 2) * fDist (fabs(omegaB(R)) / (omegaW(R) * gammaU(R)));
}
