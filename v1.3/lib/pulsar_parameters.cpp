#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <vector>

#include <sstream>
#include <algorithm>
#include <iterator>
using namespace std;

#include "constants.h"
#include "pulsar_parameters.h"

bool is_number(const string& s) {
    if (s.empty())
        return false;
    bool trigger = false;
    if (s[0] == '.') {
        return false;
    }
    for (int i = 0; i < s.size(); i ++) {
        char c = s[i];
        if (c == '.' && !trigger) {
            trigger = !trigger;
            continue;
        } else if (c == '.' && trigger) {
            return false;
        }
        if (i == 0 && c == '-') {
            continue;
        }
        if (!std::isdigit(c))
          return false;
    }
    return true;
}
double read_from_file (const string param) {
    ifstream infile("input.parameters");
    string str;
    if (infile.is_open()) {
        while (getline(infile, str)) {
            // cout << iter << ": " << str << "\n";
            istringstream iss(str);
            vector<string> words{istream_iterator<string>{iss}, istream_iterator<string>{}};
            if(words.size() < 1) continue;
            // cout << words[0] << "\n";
            if (words[0] == param) {
                if(is_number(words[1])) {
                    return atof(words[1].c_str());
                } else {
                    cout << "Error converting string to number for parameter '" << param << "' in the input file.\n";
                    exit (EXIT_FAILURE);
                }
                break;
            }
        }
    } else {
        cout << "Cannot open the file.\n";
        exit (EXIT_FAILURE);
    }
    cout << "Cannot find the given parameter '" << param << "' in the input file.\n";
    exit (EXIT_FAILURE);
}

const double B12 = read_from_file("B12"); // Surface B-field in 10^12 Gs
const double Period = read_from_file("Period"); // Rotation period in sec
const double freqGHz = read_from_file("freqGHz"); // Radiation frequency in GHz

const double lambda = read_from_file("lambda"); // Plasma multiplicity // normally 10000
const double gamma0 = read_from_file("gamma0"); // Plasma mean gamma factor // normally 50
const double f0 = read_from_file("f0"); // Polar Cap gap width // normally 0.5
const double R_em = read_from_file("R_em"); // Emission radius in star radii

const double mode = read_from_file("mode"); // 1 = O-mode & 0 = X-mode
const double fr = read_from_file("fr"); // Split-monopole parameter
const double fphi = read_from_file("fphi"); // Split-monopole parameter

const double alpha_deg = read_from_file("alpha_deg");
const double beta_deg = read_from_file("beta_deg");

const double alpha = alpha_deg * PI / 180.0; // Inclination angle in radians
const double beta = beta_deg * PI / 180.0; // Line of sight angle in radians;
const double dzeta = alpha - beta; // Minimum angle between the rotation axis and the line of sight

double PHI0;

const double B0 = B12 * 1.0e12; // Surface magnetic field in Gs
const double Omega = 2.0 * PI / Period; // Rotation frequency
const double omega = 2.0 * PI * freqGHz * 1.0e9; // Radiation circular frequency

const vector <double> vOmega{0.0, 0.0, Omega}; // Rotation vector

const double RLC = (c / Omega) / R_star;
const double RESCAPE = 1.0e3 * pow(lambda / 1.0e4, 1.0/3.0) * pow(gamma0 / 100.0, -6.0/5.0) * pow(B0 / 1.0e12, 2.0/5.0) * pow(freqGHz, -2.0/5.0) * pow(Period, -1.0/5.0);
const double ROMODE = 1.0e2 * pow(lambda / 1.0e4, 1.0/3.0) * pow(gamma0 / 100.0, 1.0/3.0) * pow(B0 / 1.0e12, 1.0/3.0) * pow(freqGHz, -2.0/3.0) * pow(Period, -1.0/3.0);
