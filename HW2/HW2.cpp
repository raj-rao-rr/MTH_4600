 /////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "Functions.h"
double YTM (double, double, int);


int main() {

    int i, j, tranches, names, name, m, maturity;
    double rho, pv0, rate, cashflow, promised_tranche_cashflow, monthtotal, temp, epsilon = 0.001;
    double **V, **L, **Y, **X, **T, **P;

    // iterating through each value of rho (correlation coefficient of 20 entities)
    // higher values of rho should create riskier cashflows, as default risks become compounded
    rho = 0.99999;

    // 20 entities in the CDO pool, allocated to 5 trenches 
    names = 20;
    tranches = 5;

    // allocating space for matrix
    T = Array (names,1);
    P = Array (tranches,1);
        
    // when Rho < 1 then we can easily perform cholesky decomposition without fault, otherwise our matrix V is not PD
    if (rho < 1.0) {
        // given that cholesky doesn't work in cases where rho = 1, since the determinant is zero we ignore the correlated normals from covariance matrix 
        temp = -50 * log(MTUniform(0));

        // since all the entities are correlated with coefficient 1, they all default at the same time  
        for (i = 1; i <= names; ++i) {
            T[i][1] = temp;
        }
    }
    else {
        // Allocate array space for V and Y.
        V = Array(names, names);
        Y = Array(names, 1);
        X = Array(names, 1);

        // Create the covariance matrix, allocating the index position accordingly
        // where: rho for (i == j) and 1 for (i != j)
        for (i = 1; i <= names; ++i) {
            for (j = 1; j <= names; ++j) {
                V[i][j] = (i == j ? 1.0 : rho);
            }
        }

        // Populate our vector X
        // Creates a vector of standard normals, provided the inverse CDF approximation, with a uniform r.v  parameter
        for (i = 1; i <= names; ++i) {
            X[i][1] = PsiInv(MTUniform(0));
        }

        // Creates our matrix L, from our PD, symmetric covariance matrix 
        // Calculate the Cholesky decomposition.  The function allocates space for L.
        L = Cholesky(V);

        // Calculates correlated gaussian random variables, with shared covariance function 
        Y = (Multiply(L, X));

        // Assign values to our default vector T
        // Note that the correlated Normals are converted into Uniforms with the Psi function
        for (i = 1; i <= names; ++i) {
            T[i][1] = -50 * log(Psi(Y[i][1]));
        }

        // Convert values of T into last month before default.
        // Taking the integer part of the default time, T*12 = m, where (int)(m) = the last month of cashflow to the CDO pool 
        for (i = 1; i <= names; ++i) {
            T[i][1] = int(T[i][1] * 12);
        }
    }
        
    // Discount Rate
    rate = 0.03;

    // Thirty year maturity in months.
    maturity = 360;

    // Promised monthly cash flow is $100.
    cashflow = 100.0;
    promised_tranche_cashflow = cashflow * names / tranches;

    // Initialize the present value of each tranche to 0.
    for (i = 1; i <= tranches; i++) {
        P[i][1] = 0.0;
    }

    // define the present value to be 0.0
    pv0 = 0.0;

    // discount the promised cash flows according to the rate and month iterated through 
    for (m = 1; m <= maturity; m ++) {
        pv0 += exp (-rate * m / 12.0) * promised_tranche_cashflow;
    }

    // Loop through the monthly payments to maturity.
    for (m = 1; m <= maturity; m ++) {

        monthtotal = 0.0;

        // calculate the cashflows per month, through maturity 
        // simultaneously check for default times based on previous calculation
        for (name = 1; name <= 20; ++name) {
            if (T[name][1] >= m){
                monthtotal += cashflow;
            }
        }

        // now we iterate through the tranches and allocate according to each tranche's promised cash flow 
        for (i = 1; i <= tranches; ++i) {
            if (monthtotal > 0) {
                // if the monthly cashflow is greater than the promised cashflow we allocate the cashflow to the respective tranche
                // we then reduce the monthly cashflow by the amount paid out to the tranche 
                if (monthtotal >= promised_tranche_cashflow) {
                    P[i][1] += exp (-rate * m / 12.0) * promised_tranche_cashflow;
                    monthtotal -= promised_tranche_cashflow;
                }
                // if the monthly cashflow is less than the promised cashflow, we allocate the remaining and fill the void with 0
                else {
                    P[i][1] += exp (-rate * m / 12.0) * monthtotal;
                    monthtotal = 0;
                }
            }
        }
    }


    // Initialize the present value of each tranche to 0 and output the values per each tranche
    for (i = 1; i <= tranches; i++) {
        printf ("For trache number %d, the present value of the expected cash flow is %8.2f.\n", i, P[i][1]);
        printf ("As a percent of the riskless PV, this is %8.2f.\n", 100.0 * P[i][1]/pv0);
        printf ("The yield-to-maturity is %8.2f percent.\n\n\n", YTM(P[i][1], promised_tranche_cashflow, maturity));
    }
    
}


////////////////////////////////////////////////////////////////////////////////
// This function calculates the yield-to-maturity (ytm) of the promised cash flow
//   of the annuity where the calculated value of this promised cash flow
//   is "pv". The ytm is calculated as a continuously compounding annual rate.
////////////////////////////////////////////////////////////////////////////////
double YTM (double pv, double cashflow, int maturity) {

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
            pv0 += cashflow * exp (-ytm * m / 12.0);
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