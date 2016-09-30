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
    alpha = alpha_deg * PI / 180.0; // Inclination angle in radians
    beta = beta_deg * PI / 180.0; // Line of sight angle in radians;
    dzeta = alpha - beta; // Minimum angle between the rotation axis and the line of sight

    vOmega[0] = 0.0;
    vOmega[1] = 0.0;
    vOmega[2] = Omega;

    cin >> R_em;

    cout << "\nR_esc: " << RESCAPE << endl;
    cout << "R_A: " << ROMODE << endl;
    cout << "R_lc: " << RLC << endl << endl;

    ostringstream temp;
    temp << (int)R_em;

    string name = "B0540+23_REM=" + temp.str();
    string path = "dats/individual/";

    string MODE;
    if (mode == 0) MODE = "X-mode";
    else MODE = "O-mode";
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

/*
    ofstream test("/Users/haykh/Dropbox/Documents/Science/Pulsars/PSR_PROG/PULSAR_PROPAGATION/v1.0/dats/test.dat");

    PHI0 = -2.0 * PI / 180.0;
    findInitPoints (PHI0);
    for (double rr = 1.0; rr <= 2000.0; rr += 1.0) {
        test << rr << " " << vBtor(rr)[0] << " " << vBtor(rr)[1] << " " << vBtor(rr)[2] << endl;
    }

    test.close();

    return 0;
*/
    //
    //    SIMULATION STARTS HERE
    //
/*
    double alphas[6] = {10.0, 25.0, 40.0, 60.0, 75.0, 80.0};
    double betas[6] = {-10.0, -7.0, -2.0, 2.0, 7.0, 10.0};

    for (int i = 0; i < 6; i ++) {
    alpha_deg = alphas[i];
    for (int j = 0; j < 6; j ++) {
    beta_deg = betas[j];
*/


    ofstream output0(path + name + "[0].dat");
    ofstream output1(path + name + "[1].dat");

    for (double phi_t = 30; phi_t >= -30; phi_t -= 0.5) { // Phase switch
        cout << "PHI: " << phi_t << endl;
        PHI0 = phi_t * PI / 180.0;
        findInitPoints (PHI0);

        double x1, x2, dep_vars[2];
        x1 = 20.0;
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

        double tau = PI * R_star / (c * omega) * integrate(dtau, 1.0, RLC);
        double gf = gFunc(pow(sin(psi_m(0.0)), 2) * RLC / NORM(vR(0.0)));
        double II0 = gf;
        double II = II0 * exp (-tau);

        double VV = II * tanh(2.0 * dep_vars[1]);
        double PA = dep_vars[0] * 180 / PI;
        output0 << phi_t << " " << II0 << " " << VV << " " << PA << endl;

        int nvar = 2, nok = 0, nbad = 0;
        double deps = 0.00001, h1 = 1.0e-14, hmin = 1.0e-15;
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
