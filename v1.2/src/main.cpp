#include <iostream>
#include <math.h>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

using namespace std;

#include "../lib/constants.h"
#include "../lib/pulsar_parameters.h"

#include "../lib/functions.h"
#include "../lib/process_functions.h"
#include "../lib/auxilary.h"

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

    string name = "test_new";
    string path = "dats/";

    struct stat st = {0};
    if (stat(path.c_str(), &st) == -1) {
        mkdir(path.c_str(), 0700);
    }

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

/*    ofstream plot(path + "betadelta.dat");
    PHI0 = 5 * PI / 180.0;
    findInitPoints(PHI0);
    for (double r = 10.0; r <= 400; r += 1.0) {
    	plot << r << " " << BetaB(r) + delta(r) << "\n";
    }
    plot.close();
    return 0;
*/
    for (double phi_t = -40; phi_t <= 40; phi_t += 1.0) { // Phase switch
        cout << "PHI: " << phi_t << endl;
        PHI0 = phi_t * PI / 180.0;
        findInitPoints (PHI0);

        double x1, x2, dep_vars[2];
        x1 = 0.0;
        // if (RESCAPE < 1000.0) x2 = 1500.0;
        x2 = 1.5 * RESCAPE;

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
        double tau = PI * R_star * integrate(dtau, x1, RLC) / (c * omega);
        double gf = gFunc(x1);
        double II0 = gf;
        double II = II0 * exp (-tau);

        double VV = II * tanh(2.0 * dep_vars[1]);
        output0 << phi_t << " " << II0 << " " << VV << " " << PA << endl;

        int nvar = 2, nok = 0, nbad = 0;
        double deps = 1.0, h1 = 1.0e-14, hmin = 1.0e-15;
        odeint(dep_vars, nvar, x1, x2, deps, h1, hmin, nok, nbad, RHS);

        VV = II * tanh(2.0 * dep_vars[1]);
        PA = dep_vars[0] * 180 / PI;

        cout << "\tI: " << II << "\n\tV: " << VV << "\n\tPA: " << -PA << endl << endl;
        output1 << phi_t << " " << II << " " << VV << " " << -PA << endl;
    }
    output0.close();
    output1.close();

	return 0;
}
