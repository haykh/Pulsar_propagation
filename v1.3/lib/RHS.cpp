#include <math.h>
#include <vector>
#include "RHS.h"
#include "constants.h"
#include "pulsar_parameters.h"
#include "process_functions.h"
#include "auxilary.h"

void RHS(double R, double *f, double *dydx) {
	double coeff = R_star * omega / (2.0 * c);

	double LL = Lambda (R);
	double QQ = Q (R);
	double BB = BetaB (R);
	double DD = delta (R);

	dydx[0] = coeff * (-LL / QQ - LL * cos(2 * f[0] - 2 * BB - 2 * DD) * sinh(2 * f[1]));
	dydx[1] = coeff * LL * sin(2 * f[0] - 2 * BB - 2 * DD) * cosh(2 * f[1]);
}
