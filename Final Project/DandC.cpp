
// Found below the main program:
void     ValueSecurity ();
void     Calibrate (double *, double);
double **AllocateLatticeArray ();
double   ValueZero (int);

#include "Functions.h"

// Global variables:
double **d;            // The state-dependent single-period discount factors.
double **V=NULL;       // State-dependent lattice values of future cash flow.

int main () {

   int n;
   double sigma, r;
   double *par, *discount, *zero, *forward;

   // Allocate memory for term structure data.
   par      = List(60);
   discount = List(60);
   zero     = List(60);
   forward  = List(60);

   // Get memory for valuation.
   V = AllocateLatticeArray ();

   // Get the user specified par curve.
   r = GetDouble ("What is the flat par interest rate in percent?... ");

   // Convert to a decimal and populate the par curve.
   r /= 100.0;
   for (n = 1; n <= 60; n++) {
      par[n] = r;
   }

   // Calculate the corresponding discount function.
   ComputeCurves (par, discount, zero, forward);

   // Get the volatility parameter.
   sigma = GetDouble ("What is the short-rate volatility in percent?... ");

   // Convert to a decimal.
   sigma /= 100.0;

   // Calibrate the Salomon binomial lattice to explain the discount function.
   Calibrate (discount, sigma);

   // Value D & C's security on this lattice.
   ValueSecurity ();

}

////////////////////////////////////////////////////////////////////////////////
// YOUR CODE GOES HERE.
////////////////////////////////////////////////////////////////////////////////
void ValueSecurity () {

    // This code should value the security and generate an output
    // file Strategy.txt for viewing the optimal withdrawal strategy
    // using WithdrawalStrategy.tex.
    
    return;

}

////////////////////////////////////////////////////////////////////////////////
// This function calibrates the state-dependent single-period discount
//    factors d[n][i] to explain the time-0 discount function discount[n].
////////////////////////////////////////////////////////////////////////////////
void Calibrate (double *discount, double sigma) {

   int n, i;
   double lower, upper, trial, gamma, v_trial, v_target, r;
   FILE *fpr;

   // Open the lattice data output files.
   fpr  = fopen ("LatticeRates.txt", "w");

   // Allocate memory for the state-dependent single-period discount factors.
   d = AllocateLatticeArray ();

   // Specify the gamma parameter with delta_t = 0.5.
   gamma = exp (sigma * sqrt(0.5));

   // Calibrate forward in time.
   for (n = 0; n < 60; n++) {

      // Bounds on rho.
      lower = 0.0;
      upper = 1.0;
      v_target = discount[n+1];

      // Iterate using interval bisection until close.
      while (1) {

         // Compute trial discount factors at period n.
         trial = (lower + upper) / 2.0;
         for (i = -n; i <= n; i += 2) {
            d[n][i] = pow (trial, pow(gamma, i));
         }

         // Value a zero maturing in period n+1 using trial rho.
         v_trial = ValueZero (n+1);

         // See if error tolerance has been met. If so, terminate search.
         if (fabs (v_trial - v_target) < 0.00000001) {
            break;
         }

         else if (v_trial > v_target) {
            // Too little discounting.
            upper = trial;
         }

         else {
            // Too much discounting.
            lower = trial;
         }

      }


      // Report the state-dependent single-period rates to a file.
      for (i = -n; i <= n; i += 2) {

         // Show only rates below 15%.
         if (d[n][i] > exp(-0.5*0.15)) {
            r = -200.0 * log(d[n][i]);
            fprintf (fpr, "%8.2f %8.2f\n", 0.5*n, r);
         }

      }

   }

   // Close the output files.
   fclose (fpr);

   printf ("View the lattice structure with Curves.tex.\n");
   Pause ();

   return;

}

////////////////////////////////////////////////////////////////////////////////
// This function values a zero coupon bond that pays $1 at period m using
//    the current d[n][i]'s for 0 <= n < m.
////////////////////////////////////////////////////////////////////////////////
double ValueZero (int m) {

   int n, i;
   double C;


   // Boundary conditions.
   for (i = -m; i <= m; i += 2) {
      V[m][i] = 0.0;
   }

   // Iterate backward recursively.
   for (n = m-1; n >= 0; n--) {
      for (i = -n; i <= n; i += 2) {

         // Calculate the cash flow.
         if (n+1 == m) {
            C = 1.0;
         } else {
            C = 0.0;
         }

         V[n][i] = d[n][i] * (0.5 * (C + V[n+1][i-1]) + 0.5 * (C + V[n+1][i+1]));

      }
   }

   return V[0][0];

}

////////////////////////////////////////////////////////////////////////////////
// Allocate memory for a lattice array x[n][i], where 0 <= n <= 60 and
//    -n <= i <= n.
////////////////////////////////////////////////////////////////////////////////
double **AllocateLatticeArray () {

   int n;
   double **x;

   x = (double **) calloc (61, sizeof (double *));
   for (n = 0; n <= 60; n++) {
      x[n] = (double *) calloc (2*n + 1, sizeof (double));
      x[n] += n;
   }

   return x;

}

