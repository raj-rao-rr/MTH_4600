#include <iostream>
#include <vector>
#include "Functions.h"


// Found below the main program:
double   ValueSecurity (int&);
void     Calibrate (double *, double);
void     ReCalibrate(double&, double&);
double **AllocateLatticeArray ();
double   ValueZero (int);
std::pair<double, double> interest_rate_risk(double&, double&, double&);


// Global variables:
double **d;            // The state-dependent single-period discount factors.
double **V=NULL;       // State-dependent lattice values of future cash flow.


int main () {

   int n, plot;
   double a, b, c, sigma, r, i, deltaR, price, function;
   std::pair<double, double> risk_calculation_100, risk_calculation_1;
   std::vector<double> price_plot, quadratic_fit_100, quadratic_fit_1;
   double *par, *discount, *zero, *forward;
   FILE *fp;

   // Allocate memory for term structure data.
   par      = List(60);
   discount = List(60);
   discount = List(60);
   zero     = List(60);
   forward  = List(60);

   // Get memory for valuation.
   V = AllocateLatticeArray ();

   ////////////////////////////////////////////////////////////////////////////////
   // Part 1
   ////////////////////////////////////////////////////////////////////////////////

   printf("Problem 1:\n");
   for (i = 0; i <= 10; ++i){ 
      // Get the user specified par curve.
      r = i; //5; //GetDouble ("What is the flat par interest rate in percent?... ");
      r /= 100.0; // Convert to a decimal and populate the par curve.

      // Get the volatility parameter.
      sigma = 25.0; //24.72958; //GetDouble ("What is the short-rate volatility in percent?... ");
      sigma /= 100.0; // Convert to a decimal.

      ReCalibrate(r, sigma);

      printf ("For r = %5.2f percent, ", r*100);
      // Value D & C's security on this lattice only when r = 5%
      if (r*100 != 5.0){
         plot = 0;
         price = ValueSecurity(plot);
         printf ("the security's value is %5.2f\n", price);
      } else {
         plot = 1;
         price = ValueSecurity(plot);
         printf ("the security's value is %5.2f\n", price);
      }

      // track all prices created for provided interest rates
      price_plot.push_back(price);
   }
   
   
   ////////////////////////////////////////////////////////////////////////////////
   // Part 2 & 3 (Duration & Convexity)
   ////////////////////////////////////////////////////////////////////////////////
   
   printf("\n\nProblem 2:\n");

   r = 5.0;             // in percent
   r /= 100.0;          // nominal

   sigma = 25.0;
   sigma /= 100.0;

   // the 100bps shift in interest rates
   deltaR = 100.0;    // in bps
   deltaR /= 10000.0; // nominal

   risk_calculation_100 = interest_rate_risk(r, deltaR, sigma);
   
   // the 1bps shift in interest rates
   deltaR = 1.0;      // in bps
   deltaR /= 10000.0; // nominal

   risk_calculation_1 = interest_rate_risk(r, deltaR, sigma);

   printf("For r = %5.2f percent ", r * 100.0);
   printf("and deltaR = %5.2f bps, ", 100.0);
   printf("the calculated duration is %5.2f\n", risk_calculation_100.first);

   printf("For r = %5.2f percent ", r * 100.0);
   printf("and deltaR = %5.2f bps, ", 1.0);
   printf("the calculated duration is %5.2f\n", risk_calculation_1.first);

   printf("\n\nProblem 3:\n");

   printf("For r = %5.2f percent ", r * 100.0);
   printf("and deltaR = %5.2f bps, ", 100.0);
   printf("the calculated convexity is %5.2f\n", risk_calculation_100.second);
  
   printf("For r = %5.2f percent ", r * 100.0);
   printf("and deltaR = %5.2f bps, ", 1.0);
   printf("the calculated convexity is %5.2f\n", risk_calculation_1.second);


   ////////////////////////////////////////////////////////////////////////////////
   // Part 4
   ////////////////////////////////////////////////////////////////////////////////
  
   printf("\n\nProblem 4:\n");
   printf("Price at 5.0 is %8.4f\n", price_plot[5]);
   // quadratic fit for 100 bps estimates
   a = 0.5 * risk_calculation_100.second * price_plot[5];       // half the convexity value at 5% rate (100bps)
   b = risk_calculation_100.first * -price_plot[5];             // duration value for bond at 5% rate (100bps)
   c = price_plot[5];                                           // fair value price at 5% rate

   for (i = 0; i <= 10; ++i) {
       function = a * (i - 5) * (i - 5) + b * (i - 5) + c;
       quadratic_fit_100.push_back(function);
   }
    
   // quadratic fit for 1 bps estimate 
   a = 0.5 * risk_calculation_1.second * price_plot[5];         // half the convexity value at 5% rate (1bps)
   b = risk_calculation_1.first * -price_plot[5];               // duration value for bond at 5% rate (1bps)
   c = price_plot[5];                                           // fair value price at 5% rate

   for (i = 0; i <= 10; ++i) {
       function = a * (i - 5) * (i - 5) + b * (i - 5) + c;
       quadratic_fit_1.push_back(function);
   }

   // Output the vector for each of the quadratic estimators 
   printf("Quadratic fit vector (100bps)\n");
   for(int i=0; i <= 10; i++){
      std::cout << quadratic_fit_100[i] << ' ';
   }
   printf("\n\nQuadratic fit vector (1bps)\n");
   for (int i = 0; i <= 10; i++) {
       std::cout << quadratic_fit_1[i] << ' ';
   }

   // prices plot vector being used 
   fp = fopen ("prices.txt", "w");

   for(int i=0; i < price_plot.size(); i++){
      fprintf (fp, "\\put {$\\scriptscriptstyle\\bullet$} at %d  %6.2f\n", 10-i, price_plot.at(i));
   }
   fclose (fp);

   // quadratic fit for 1 bps measure
   fp = fopen ("quadraticfit_1.txt", "w");
   for(int i=0; i < quadratic_fit_1.size(); i++){
      fprintf (fp, "%d  %6.2f\n", i, quadratic_fit_1.at(i));
   }
   fclose (fp);

   // quadratic fit for 100 bps measure 
   fp = fopen ("quadraticfit_100.txt", "w");
   for(int i=0; i < quadratic_fit_1.size(); i++){
      fprintf (fp, "%d  %6.2f\n", i, quadratic_fit_100.at(i));
   }
   fclose (fp);

   printf("\n");
   Exit();
}

// computes both the duration and convexity of a bond
std::pair<double, double> interest_rate_risk(double& r, double& deltaR, double& sigma) {
    std::pair<double, double> dur_cov;
    double new_r, pBase, pPlus, pMinus, dur, cov;
    int plot = 0;

    // pBase (the base line formula for the computation)
    ReCalibrate(r, sigma);
    pBase = ValueSecurity(plot);

    // pPlus (adding a factor of delta r)
    new_r = r + deltaR;
    ReCalibrate(new_r, sigma);
    pPlus = ValueSecurity(plot);

    // pMinus (subtracting a factor of delta r)
    new_r = r - deltaR;
    ReCalibrate(new_r, sigma);
    pMinus = ValueSecurity(plot);

    // compute both the duration and convexity
    dur = (-1.0 / pBase) * (pPlus - pMinus) / (2.0 * deltaR);
    cov = (1.0 / pBase) * (pPlus - 2.0 * pBase + pMinus) / (deltaR * deltaR);

    // assigning the pair (tuple) to values for the duration and convexity
    dur_cov.first = dur; dur_cov.second = cov;

    return dur_cov;
}


void ReCalibrate(double& r, double& sigma){
   double *par, *discount, *zero, *forward;

   par      = List(60);
   discount = List(60);
   zero     = List(60);
   forward  = List(60);

   for (int n = 1; n <= 60; n++) {
      par[n] = r;
   }

   ComputeCurves (par, discount, zero, forward);
   Calibrate (discount, sigma);
}


////////////////////////////////////////////////////////////////////////////////
// YOUR CODE GOES HERE.
////////////////////////////////////////////////////////////////////////////////
double ValueSecurity (int& plot) {

   int n, i, m, called;
   double r, accumulated_cash;
   FILE *fp;

   if (plot == 1){
      // Open output file to contain the call option exercise strategy.
      fp = fopen ("Strategy.txt", "w");
   }

   // Get the bond's maturity in years and convert to semiannual periods.
   m  = 30; //GetInteger ("What is the bond's maturity in years (integer <= 30, please)?... ");
   m *= 2;

   // Value $100 face amount.
   // Boundary conditions.
   for (i = -m; i <= m; i += 2) {
      //         60 * (300 - 5 * (60-1) / 2) = 9150
      V[m][i] = (m) * (300 - 5.0 * (m-1) / 2.0);
   }

   // Iterate backward recursively.
   for (n = m-1; n >= 0; n--) {

      accumulated_cash = (n+1) * (300 - (5.0 * n / 2.0) );

      for (i = -n; i <= n; i += 2) {

         // Value today if the security remains outstanding until period n+1.
         V[n][i] = d[n][i] * (V[n+1][i-1] + V[n+1][i+1]) / 2.0 ;

         //*********************************************************************
         // See if the security should be called.
         // If the value of the accumulated cash flow is higher than
         // the value today if the security remains outstanding until period n+1,
         // the security should be called.
         if (accumulated_cash > V[n][i]){
            V[n][i] = accumulated_cash;
            called = 1;
         } else {
            called = 0;
         }
         //*********************************************************************

         // Report the option exercise strategy for viewing with TeX software.
         // Show only for short rates r <= 15%, which corresponds to
         //    d[n][i] >= exp(-0.15*0.5) = 0.927743.
         if (plot == 1){
            if (d[n][i] >= 0.927743) {
               r = -200.0 * log (d[n][i]);
               if (called) {
                  fprintf (fp, "\\put{$\\scriptscriptstyle\\bullet$} at %5.1f  %6.2f\n", 0.5*n, r);
               } else {
                  fprintf (fp, "\\put{$\\scriptscriptstyle\\circ$} at %5.1f  %6.2f\n", 0.5*n, r);
               }
            }
         }

      }
   }

   if (plot == 1){
      // Close the output file.
      fclose (fp);
   }
   
   return V[0][0];
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

   // printf ("View the lattice structure with Curves.tex.\n");
   //Pause ();

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

