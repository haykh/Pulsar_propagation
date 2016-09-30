#include <math.h>
#include <vector>
using namespace std;

#include "constants.h"
#include "pulsar_parameters.h"

const double B12 = 1.0;
const double Period = 1.0; // Rotation period in sec // normally 1.0
const double freqGHz = 1.0; // Radiation frequency in GHz // normally 1.0

const double lambda = 500.0; // Plasma multiplicity // normally 10000
const double gamma0 = 25.0; // Plasma mean gamma factor // normally 50
const double f0 = 0.6; // Polar Cap gap width // normally 0.5
const double R_em = 50.0; // Emission radius in star radii

const double mode = 1; // 1 = O-mode & 0 = X-mode
const double fr = 1.0; // Split-monopole parameter
const double fphi = 1.0; // Split-monopole parameter
const double BMULT = 0.0; // Toroidal field multiplier

const double alpha_deg = 70.0;
const double beta_deg = 7.0;

double alpha; // Inclination angle in radians
double beta; // Line of sight angle in radians;
double dzeta; // Minimum angle between the rotation axis and the line of sight

double PHI0;

vector <double> vOmega(3); // Rotation vector

const double B0 = B12 * 1.0e12; // Surface magnetic field in Gs // normally 10^12
const double Omega = 2.0 * PI / Period; // Rotation frequency
const double omega = 2.0 * PI * freqGHz * 1.0e9; // Radiation circular frequency

const double RLC = (c / Omega) / R_star;
const double RESCAPE = 1.0e3 * pow(lambda / 1.0e4, 1.0/3.0) * pow(gamma0 / 100.0, -6.0/5.0) * pow(B0 / 1.0e12, 2.0/5.0) * pow(freqGHz, -2.0/5.0) * pow(Period, -1.0/5.0);
const double ROMODE = 1.0e2 * pow(lambda / 1.0e4, 1.0/3.0) * pow(gamma0 / 100.0, 1.0/3.0) * pow(B0 / 1.0e12, 1.0/3.0) * pow(freqGHz, -2.0/3.0) * pow(Period, -1.0/3.0);
