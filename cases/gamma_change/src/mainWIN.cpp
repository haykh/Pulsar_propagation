#include <iostream>
#include <math.h>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
using namespace std;

#include "../lib/constants.h"
#include "../lib/pulsar_parameters.h"

#include "../lib/functions.h"
#include "../lib/geometric_functions.h"
#include "../lib/physical_functions.h"

#include "../lib/integrator.h"
#include "../lib/RHS.h"

/*
    Default libraries
*/
#include "../lib/diffeqsolver.h"

void displayVector (vector <double> a) {
    cout << endl << a[0] << endl << a[1] << endl << a[2] << endl;
}

int main() {
    cout << "gamma = ";
    cin >> gamma0;

    RESCAPE = 1.0e3 * pow(lambda / 1.0e4, 1.0/3.0) * pow(gamma0 / 100.0, -6.0/5.0) * pow(B0 / 1.0e12, 2.0/5.0) * pow(freqGHz, -2.0/5.0) * pow(Period, -1.0/5.0);
    ROMODE = 1.0e2 * pow(lambda / 1.0e4, 1.0/3.0) * pow(gamma0 / 100.0, 1.0/3.0) * pow(B0 / 1.0e12, 1.0/3.0) * pow(freqGHz, -2.0/3.0) * pow(Period, -1.0/3.0);

    cout << "\nR_esc: " << RESCAPE << endl;
    cout << "R_A: " << ROMODE << endl;
    cout << "R_lc: " << RLC << endl << endl;

    alpha = alpha_deg * PI / 180.0; // Inclination angle in radians
    beta = beta_deg * PI / 180.0; // Line of sight angle in radians;
    dzeta = alpha - beta; // Minimum angle between the rotation axis and the line of sight

    vOmega[0] = 0.0;
    vOmega[1] = 0.0;
    vOmega[2] = Omega;

    string MODE;
    if (mode == 0) MODE = "X-mode";
    else MODE = "O-mode";

    ostringstream temp;
    temp << (int)gamma0;

    string name = "gamma0=" + temp.str();
    string path = "dats/";

    ofstream outputData(path + name + ".dat");
    outputData
        << "alpha = " << alpha_deg
        << "\nbeta = " << beta_deg
        << "\n\nPeriod = " << Period
        << "\nB12 = " << B12
        << "\nfGHz = " << freqGHz
        << "\n\nlambda = " << lambda
        << "\ngamma0 = " << gamma0
        << "\nf0 = " << f0
        << "\nR_em = " << R_em
        << "\n\n" + MODE
        << "\nfr = " << fr
        << "\nfphi = " << fphi
        << "\nBtor = " << BMULT
        << "\n\n\nR_LC = " << RLC
        << "\nR_escape = " << RESCAPE
        << "\nR_A = " << ROMODE;
    outputData.close();

//
//    SIMULATION STARTS HERE
//

    ofstream output0(path + name + "[0].dat");
    ofstream output1(path + name + "[1].dat");

    for (double phi_t = 20; phi_t >= -20; phi_t -= 0.5) { // Phase switch
        cout << "PHI: " << phi_t << endl;
        PHI0 = phi_t * PI / 180.0;
        findInitPoints (PHI0);

        double x1, x2, dep_vars[2];
        x1 = 10;
        if (RESCAPE < 1000.0) x2 = 1000.0;
        else x2 = RESCAPE;

        /*Initial values*/
        if (mode == 0) { // X-mode
            dep_vars[0] = BetaB(x1) + delta(x1) + PI / 2.0;
            dep_vars[1] = Arcsinh(1.0 / Q(x1)) / 2.0;
        } else { // O-mode
            dep_vars[0] = BetaB(x1) + delta(x1);
            dep_vars[1] = Arcsinh(-1.0 / Q(x1)) / 2.0;
        }
        /*--------------*/

        double PA = dep_vars[0] * 180 / PI;
        double tau = PI * R_star * integrate(dtau, 1.0, RLC) / (c * omega);
        double gf = gFunc(pow(sin(psi_m(0.0)), 2) * RLC / NORM(vR(0.0)));
        double II0 = gf;
        double II = II0 * exp (-tau);

        double VV = II * tanh(2.0 * dep_vars[1]);
        output0 << phi_t << " " << II0 << " " << VV << " " << PA << endl;

        int nvar = 2, nok = 0, nbad = 0;
        double deps = 1.0, h1 = 1.0e-14, hmin = 1.0e-15;
        odeint(dep_vars, nvar, x1, x2, deps, h1, hmin, nok, nbad, RHS);

        VV = II * tanh(2.0 * dep_vars[1]);
        PA = dep_vars[0] * 180 / PI;

        cout << "\tI: " << II << "\n\tV: " << VV << "\n\tPA: " << PA << endl << endl;
        output1 << phi_t << " " << II << " " << VV << " " << PA << endl;
    }
    output0.close();
    output1.close();

	return 0;
}
