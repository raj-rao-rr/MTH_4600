////////////////////////////////////////////////////////////////////////////////
// Starter code for Homework 6.
////////////////////////////////////////////////////////////////////////////////
#define _USE_MATH_DEFINES

// These functions are found below.
void GetData ();
void Report ();

// These global variables are used to reduce function arguments.
double *R, *R2, *sigma;

#include "Functions.h"
#include <iostream>
#include <cmath>


int main () {

   double *s, *p, *ptilde, *g, *annual_vol;
   double sigma0, alpha, sqrt2pi, k_hat_max, Z, max_probT, current_probT;
   const double pi = M_PI;
   double *AllocateMemory();
   int index, t, k, i;

   // Allocate necessary memory.
   s = AllocateMemory ();
   p = AllocateMemory ();
   ptilde = AllocateMemory ();
   g = AllocateMemory ();
   annual_vol = (double*) calloc(1260, sizeof(double));

   // Read in the time series data.
   GetData ();

   // Assign the free parameters 
   sigma0 = 0.25;
   alpha = 0.03;

   // Compute the value of sqrt(2*pi) once and for all
   sqrt2pi = sqrt(2 * pi);

   // P_D_0 [E] = P [E | Omega] = P [E]
   p[0] = 1;

   // Increment the time 
   for (t = 1; t <= 1258; ++t) {
        
       // Initialize certain variables
       k_hat_max = 0.0;   // k value that maximizes the conditional probability 
       max_probT = 0.0;   // Maximum value for the conditional distribution at time t  
       Z = 0.0;           // Sum of  g_k through all states k with positive probability

       // Cycle through states k with positive probability
       for (k = -t; k <= t; k += 2) {

           // Index for all the positions of k
           index = k;

           //////////////////////////////
           // 1. Volatility Estimation 
           //////////////////////////////

           // Estimate the conditional distribution of the location of the walk at time t from D_{t-1}
           // p{k} = (p{k-1}+p{k+1}) / 2
           current_probT = 0.5 * (p[index + 1] + p[index - 1]);

           // Determine which value of k  maximizes the probability P_D_{t-1} [W_t = k] 
           if (current_probT > max_probT) {
               k_hat_max = k;
               max_probT = current_probT;
           }

           //////////////////////////////
           // 2. Bayesian Updating 
           //////////////////////////////

           // Calculate sk
           s[index] = sigma0 * exp(alpha * k);
           // Calculate gk
           g[index] = 1 / (sqrt2pi * s[index]) * exp(-R2[t] / (2 * s[index] * s[index])) * current_probT;
           // Update Z, the sum of g across all states k
           Z += g[index];

       }

       // The daily volatility at k-max (this is our Volatility estimate) 
       sigma[t] = sigma0 * exp(alpha * k_hat_max);

       // Annualized volatility used to create a timeseries plot of the volatility
       annual_vol[t] = sqrt(252) * sigma[t]; 
       printf("At time %4d, k_hat is %3.0f and the annualized volatility is -> %3.4f\n", t, k_hat_max, sqrt(252)*sigma[t]);

       // Reassign the values of p with new posterior p-tilde
       for (i = -t; i <= t; i += 2) {
           // ptilde as constructed from the conditional probability P_D_{t-1} [R_t = r_t | W_t = k] * P_D_{t-1} [W_t = k] / P_D_{t-1} [R_t = r_t]
           ptilde[i] = g[i] / Z; 
           // Current ptilde's become p's at t+1
           p[i] = ptilde[i];
       }

   }

   // Create TeX files for viewing results.
   Report ();

   Exit();
}


////////////////////////////////////////////////////////////////////////////////
// Allocate memory for an array with indices from -1260 to +1260.
// This is a little more than needed.
////////////////////////////////////////////////////////////////////////////////

double *AllocateMemory () {

   double *x;

   x = (double *) calloc (2*1260 + 1, sizeof (double));

   x += 1260;

   return x;

}

////////////////////////////////////////////////////////////////////////////////
// Read in a daily time series of stock price returns R[t] 1 <= t <= 1258.
////////////////////////////////////////////////////////////////////////////////

void GetData () {

   int t;
   char input[100];
   FILE *fp;

   fp = fopen ("/Users/Paul/Desktop/MTH 4600/Visual Studio Code Folder/HW6 - new/XOM5YrsDaily.txt", "r");

   // Read in the file description.
   fgets (input, 99, fp);

   // Allocate memory for the data, initialized to 0.
   R     = List (1258);
   R2    = List (1258);
   sigma = List (1258);

   // Now read in the data.
   for (t = 1; t <= 1258; t++) {
      fgets (input, 99, fp);
      sscanf (input, "%lf", R+t);
      R2[t] = R[t] * R[t];
   }

   fclose (fp);

   return;

}

////////////////////////////////////////////////////////////////////////////////
// Generate some output files.
////////////////////////////////////////////////////////////////////////////////

void Report () {

   int t;
   double annualizedVol, Z;
   FILE *fp1, *fp2;


   fp1 = fopen ("/Users/Paul/Desktop/MTH 4600/Visual Studio Code Folder/HW6 - new/HistEmpVols.txt", "w");
   fp2 = fopen ("/Users/Paul/Desktop/MTH 4600/Visual Studio Code Folder/HW6 - new/StandardizedXOM.txt", "w");

   // Start at day 51.
   for (t = 51; t <= 1258; t++) {

      annualizedVol = sqrt (252.0) * sigma[t];
      fprintf (fp1, "%4d %10.2f\n", t, annualizedVol);

      // Standardize the data and add to a histogram.
      // Z = PsiInv(MTUniform(0));               // **** Get rid of this line,
      Z = R[t] / sigma[t];                  // and un-comment this line.
      NormalHistogram (Z, 20, 0);
      fprintf (fp2, "%4d %10.2f\n", t, Z);

   }

   fclose (fp1);
   fclose (fp2);

   // Create histogram TeX file.
   NormalHistogram (0, 20, 1);

   return;

}






