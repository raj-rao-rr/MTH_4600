////////////////////////////////////////////////////////////////////////////////
// Paul HW2
// 

#include "Functions.h"
double YTM (double, double, int);

int main() {

    int i, j, tranches, names, name, m, maturity;
    double rho, pv0, rate, cashflow, promised_tranche_cashflow, monthtotal;
    double **V, **L, **Y, **X, **T, **P;

    for (rho=0.0; rho < 1.0; rho = rho + 0.1){
        
        names = 20;
        tranches = 5;

        //printf ("I am computing the Cholesky decomposition for the matrix V\n");

        // Allocate array space for V and Y.
        V = Array (names,names);
        Y = Array (names,1);
        X = Array (names,1);
        T = Array (names,1);
        P = Array (tranches,1);

        // Assign values to the V[i][j].
        for (i = 1; i <= names; i++) {
            for (j = 1; j <= names; j++) {
                V[i][j] = (i == j ? 1.0 : rho);
            }
        }

        // Assign values to the X[i][1].
        for (i = 1; i <= names; i++) {
            X[i][1] = MTUniform(0);
        }
            
        // Show the matrix V.
        //printf ("V:\n");
        //Show (V);

        // Show the matrix X.
        //printf ("X:\n");
        //Show (X);

        // Calculate the Cholesky decomposition.  The function allocates space for L.
        L = Cholesky (V);

        // Show the matrix L.
        //printf ("L:\n");
        //Show (L);

        // Show the product L times L transpose.
        //printf ("L * (L transpose):\n");
        //Show (Multiply (L, Transpose(L)));

        // Show the product Y times L.
        Y = (Multiply (L, X));
        //printf ("L * X:\n");
        //Show (Multiply (L, X));
        //Show (Y);

        // Assign values to the T[i][1].
        // Note that the correlated Normals are converted into Uniforms with the Psi function
        for (i = 1; i <= names; i++) {
            T[i][1] = -50 * log( Psi(Y[i][1]) );
        }

        // Show the matrix T.
        //printf ("T:\n");
        //Show (T);

        // Convert values of T into last month before default.
        for (i = 1; i <= names; i++) {
            T[i][1] = int( T[i][1] * 12 );
        }

        // Show the matrix T.
        //printf ("T:\n");
        //Show (T);


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

        // Show the matrix P.
        //printf ("P:\n");
        //Show (P);

        pv0 = 0.0;

        for (m = 1; m <= maturity; m ++) {
            pv0 += exp (-rate * m / 12.0) * promised_tranche_cashflow;
        }

        // Loop through the monthly payments to maturity.
        for (m = 1; m <= maturity; m ++) {

            monthtotal = 0.0;

            for (name = 1; name <= 20; ++name) {
                if (T[name][1] >= m){
                    monthtotal += cashflow;
                }
            }

            for (i = 1; i <= tranches; ++i) {
                if (monthtotal > 0) {
                    if (monthtotal >= promised_tranche_cashflow) {
                        P[i][1] += exp (-rate * m / 12.0) * promised_tranche_cashflow;
                        monthtotal -= promised_tranche_cashflow;
                    }
                    else {
                        P[i][1] += exp (-rate * m / 12.0) * monthtotal;
                        monthtotal = 0;
                    }
                }
            }
        }

        // Show the matrix P.
        //printf ("P:\n");
        //Show (P);


        // Initialize the present value of each tranche to 0.
        printf ("For Rho = %8.2f:\n",rho);
        for (i = 1; i <= tranches; i++) {
            printf ("For trache number %d, the present value of the expected cash flow is %8.2f.\n", i, P[i][1]);
            printf ("As a percent of the riskless PV, this is %8.2f.\n", 100.0 * P[i][1]/pv0);
            printf ("The yield-to-maturity is %8.2f percent.\n\n\n", YTM(P[i][1], promised_tranche_cashflow, maturity));
        }
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