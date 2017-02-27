#include <vector>
#include <math.h>
#include "functions.h"
using namespace std;

double SCALAR (vector <double> vec1, vector <double> vec2) {
    return vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
} // Scalar product
vector <double> CROSS (vector <double> vec1, vector <double> vec2) {
    vector <double> temp(3);
    temp[0] = -vec1[2] * vec2[1] + vec1[1] * vec2[2];
    temp[1] = vec1[2] * vec2[0] - vec2[2] * vec1[0];
    temp[2] = -vec1[1] * vec2[0] + vec2[1] * vec1[0];
    return temp;
} // Vector product
double NORM (vector <double> vec) {
    return sqrt(SCALAR(vec, vec));
} // Length of vector
vector <double> SUM (vector <double> vec1, vector <double> vec2) {
    vector <double> temp(3);
    temp[0] = vec1[0] + vec2[0];
    temp[1] = vec1[1] + vec2[1];
    temp[2] = vec1[2] + vec2[2];
    return temp;
} // Sum

vector <double> TIMES (double a, vector <double> vec) {
    vector <double> temp(3);
    temp[0] = vec[0] * a;
    temp[1] = vec[1] * a;
    temp[2] = vec[2] * a;
    return temp;
} // Multiply vector by double
vector <double> NORMALIZE (vector <double> vec) {
    return TIMES(1.0 / NORM(vec), vec);
} // Unitize given vector

double ANGLE (vector <double> vec1, vector <double> vec2) {
    return acos(SCALAR(NORMALIZE(vec1), NORMALIZE(vec2)));
} // Angle
double Arcsinh(double x) {
    return log(x + sqrt(pow(x, 2.0) + 1.0));
}