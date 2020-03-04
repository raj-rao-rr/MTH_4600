
double YTM (double, double, int);

#include "Functions.h"
////////////////////////////////////////////////////////////////////////////////
// MTH 4600 HW2
// Created by Group of Atreish, Ben, Cheng, Karen
// Last edited: 3/1/2020 10:34 pm
////////////////////////////////////////////////////////////////////////////////

int main() {

	double *U = (double *)calloc(21, sizeof(double)); //Independent Uniform(0,1)
	double *X = (double *)calloc(21, sizeof(double)); //Independent Standard Normal
	double *Y = (double *)calloc(21, sizeof(double)); //Dependent Standard Normal
	double *U_new = (double *)calloc(21, sizeof(double)); //Dependent Uniform(0,1)
	double *T = (double *)calloc(21, sizeof(double)); //Dependent Time of Default

	double **V = Array(20, 20); // Covariance Matrix
	double **L = Array(20, 20); // Lower Triangular Matrix from Cholesky Decomposition

	double totalCashflowPerMonth;
	double *PV = (double *)calloc(6, sizeof(double));

	double riskfreeRate = 0.03;
	double discountFactor;

	int numberOfIteration = 50000;
	double rho = 0.0;

// Loop through all values of rho from 0 to 1 with step of 0.1
while (rho < 1.000001) {

	// Initialization
	for (int i = 1; i < 21; i++) {
		U[i] = 0;
		X[i] = 0;
		Y[i] = 0;
		U_new[i] = 0;
		T[i] = 0;
	}

	for (int i = 1; i < 6; i++) {
		PV[i] = 0;
	}

	// Construction of Covariance Matrix V
		if (rho < 0.999999) {
			for (int i = 1; i < 21; i++) {
				for (int j = 1; j < 21; j++) {
					if (i == j) {
						V[i][j] = 1;
					}
					else {
						V[i][j] = rho;
					}
				}
			}

			// Compute the Cholesky decomposition V = LL^T for the symmetric PD matrix V.
			L = Cholesky(V);
		}

		//Monte Carlo simulation for EPV for 5 tranches
		for (int n = 1; n <= numberOfIteration; n++) {

			// Generate Random Variables
			if (rho < 0.999999) {
				// Independent uniform(0,1) and standard normals.
				for (int i = 1; i < 21; i++) {
					U[i] = MTUniform(0);
					X[i] = PsiInv(U[i]);
				}

				// Dependent Standard Normal Y = LX
				for (int i = 1; i < 21; i++) {
					Y[i] = 0;
					for (int k = 1; k < 21; k++) {
						Y[i] += L[i][k] * X[k];
					}
				}

				// Dependent Uniform(0,1) and Dependent Time of Default
				for (int i = 1; i < 21; i++) {
					U_new[i] = Psi(Y[i]);
					T[i] = -50 * log(U_new[i]);
				}

			}
			else { // Seperate Case when rho = 1 (Cholesky doesn't work on this case!)
				double U = MTUniform(0);
				T[1] = -50 * log(U);
				for (int i = 2; i < 21; i++) {
					T[i] = T[1];
				}
			}

			// Generate total present values of 50000 simulation for each tranche
			for (int m = 1; m < 361; m++) {
				//int m = 360;

				totalCashflowPerMonth = 0.0; //Zero out
				discountFactor = exp(-riskfreeRate * m / 12.0);

				// Monthly payments from 20 entities
				for (int k = 1; k < 21; k++) {
					if (T[k] > m / 12.0 ) {
						totalCashflowPerMonth += 100.0;
					}
				}

				// Monthly Paid-out to each tranche
				int i = 1;
				while (totalCashflowPerMonth >= 400.0) {
					PV[i] += 400.0 * discountFactor;
					totalCashflowPerMonth = totalCashflowPerMonth - 400.0;
					i++;
				}
				PV[i] += totalCashflowPerMonth * discountFactor;
			}

		}

		// Printing results
		double EPV;
		double ytm;
		char* creditRating;

		printf("For rho = %f", rho);
		printf(": \n");

		for (int i = 1; i < 6; i++) {
			EPV = PV[i] / numberOfIteration; //compute average PV over the 50000 simulations
			ytm = YTM(EPV, 400.0, 360);

			// Assign Credit Rating
			if (ytm <= 4.0) {
				creditRating = "AAA";
			}
			else {
				if (ytm <= 4.5) {
					creditRating = "AA";
				}
				else {
					if (ytm <= 5.0) {
						creditRating = "A";
					}
					else {
						if (ytm <= 6.0) {
							creditRating = "BBB";
						}
						else {
							creditRating = "Speculative";
						}
					}
				}
			}

			printf("    Tranche %d has EPV of %f", i, EPV);
			printf(", YTM of %8.2f",ytm );
			printf("%% and credit rating of %s \n", creditRating);
		}
		printf("-----------------------------");
		printf("\n");

		rho += 0.1;
	}
   Pause ();

}

////////////////////////////////////////////////////////////////////////////////
// This function calculates the yield-to-maturity (ytm) of the promised cash flow
//   of the annuity where the calculated value of this promised cash flow
//   is "pv". The ytm is calculated as a continuously compounding annual rate.
////////////////////////////////////////////////////////////////////////////////
double YTM(double pv, double cashflow, int maturity) {

	int m;
	double pv0, ytm = 0.0, step = 0.01;

	// Calculate the ytm to within 0.000001, i.e., 0.0001% or 0.01 bps
	while (step > 0.000001) {

		// Keep increasing the trial ytm until the resulting present value (pv0)
		//  is too small (less than pv).
		while (1) {

			// Calculate the present value with a discount rate of "ytm".
			pv0 = 0;
			for (m = 1; m <= maturity; m++) {
				pv0 += cashflow * exp(-ytm * m / 12.0);
			}

			// Is ytm now too big? If not, increase it further by the amount "step".
			if (pv0 > pv) {
				ytm += step;
			}

			// If so, reduce it by "step" and break out of the "while (1)" loop.
			// Reduce the step size for the next iteration.
			else {
				ytm -= step;
				step /= 10.0;
				break;
			}

		}

	}

	// Return the yield-to-maturity as a percent.
	return 100.0 * ytm;

}


