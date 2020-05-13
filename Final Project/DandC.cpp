
// Found below the main program:
void     ValueSecurity (double *);
void     Calibrate (double *, double);
double **AllocateLatticeArray ();
double   ValueZero (int);
double duration(double&, double&, double&, double&);
double convexity(double&);

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
   ValueSecurity (par);

}


////////////////////////////////////////////////////////////////////////////////
// Calculating the fixed income "greeks" 
////////////////////////////////////////////////////////////////////////////////
double duration(double& pBase, double& pPlus, double& pMinus, double& deltaR) {
    double dur;
    dur = (-1 / pBase) * (pPlus - pMinus) / (2 * deltaR);
    return dur; 
}


double convexity(double& pBase, double& pPlus, double& pMinus, double& deltaR) {
    double conv;
    conv = (1 / pBase) * (pPlus - 2*pBase +pMinus) / (deltaR * deltaR);
    return conv;
}

////////////////////////////////////////////////////////////////////////////////
// Valuing the price of the security
////////////////////////////////////////////////////////////////////////////////
void ValueSecurity (double* par) {

    // This code should value the security and generate an output
    // file Strategy.txt for viewing the optimal withdrawal strategy
    // using WithdrawalStrategy.tex.
    int m, t, i, called;
    double r, C = 0.0; 
    FILE* fps;

    // Open the optimal withdrawal strategy files.
    fps = fopen("Strategy.txt", "w");

    // Get the bond's maturity in years and convert to semiannual periods.
    m = 2 * 30;

    // Report the par rate for that maturity. (*NOTE WE HAVE A FLAT PAR YEILD CURVE*)
    printf("The par yield for that maturity is %6.3f\n", 100.0 * par[m]);

    // Boundary conditions.
    for (i = -m; i <= m; i += 2) {
        V[m][i] = 0.0;
    }

    // Iterate backward recursively in time  
    for (t = m; t >= 0; --t) {

        // checks all nodes in pairs at time t
        for (i = -t; i <= t; i += 2) {

            // Discount the cash flow at time t
            C = 300 - 5 * (m);

            // Value today if the bond remains outstanding until period n+1.
            V[t][i] = d[t][i] * (0.5 * (C + V[t + 1][i - 1]) + 0.5 * (C + V[t + 1][i + 1]));

            //////////////////////////////////
            // Excerising fixed income contract   
            //////////////////////////////////
            if (V[t][i] > 100.0) {
                V[t][i] = 100.0;
                called = 1;
            }
            else {
                called = 0;
            }
            
            //////////////////////////////////
            // Report the option withdraw strategy for viewing with TeX software.
            // Show only for short rates r == 5%, which corresponds to
            //    d[n][i] == exp(-0.05*0.5)
            //////////////////////////////////
            if (d[t][i] == exp(-0.05*0.5)) {
                r = -200.0 * log(d[t][i]);
                if (called) {
                    fprintf(fps, "\\put{$\\scriptscriptstyle\\bullet$} at %5.1f  %6.2f\n", 0.5 * t, r);
                }
                else {
                    fprintf(fps, "\\put{$\\scriptscriptstyle\\circ$} at %5.1f  %6.2f\n", 0.5 * t, r);
                }
            }

        }
    }

    // Close the output file.
    fclose(fps);

    printf("\n");
    printf("The bond's value is %5.2f\n", V[0][0]);
    Exit();

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

