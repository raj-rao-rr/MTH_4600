////////////////////////////////////////////////////////////////////////////////
// This is starter code for Project 5.
////////////////////////////////////////////////////////////////////////////////

// This function is found below.
void GetData ();

// "ticker" is a global variable.
char **ticker;

// The covariance matrix is a global variable.
double **V;

#include "Functions.h"

// Calculates the variance of the given portfolio, provided weights (wt) and covariance matrix (V)
//double* Variance(double** &V, double** &wt) {
//    return  Multiply(Multiply(wt, V), Transpose(wt));
//}


int main () {

   int i, seed;

   // Read in and display the tickers and covariance matrix.
   GetData ();

   // Show the 50 tickers.
   for (i = 1; i <= 50; i++) {
      printf (ticker[i]);
   }
   Pause ();

   // Show the covariance matrix.
   Show (V);

   // Seed the RNG.
   seed = GetInteger ("RNG Seed?... ");
   MTUniform (seed);

   /////////////////////////////////////////////////////////////////////////////
   // Your Metropolis algorithm starts here...

   // problem 1. 
   double **I, **actual_mv, **c;
   double** e = Array(1, 50);

   // define the e-array equal to [1, 1, ... , 1]^T
   for (int i = 0; i < 50; ++i) {
       e[1][i] = 1;
   }

   // compute the inverted matrix V^{-1}
   I = Invert(V);

   // compute the value of c = e^T*(V^{-1})*e
   c = Multiply(Multiply(e,I), Transpose(e));

   // the actual minimum variance portfolio for the provided covariance matrix 
   actual_mv = Multiply(ScalarMultiple((1/c[1][1]), I), Transpose(e));

   printf("The actual minimum variance portoflio\n");
   for (int i = 1; i < 50; ++i) {
       printf("%8.4f\n", actual_mv[i][1]);
   }

   // variance and weight array 
   double var;
   double** wt = Array(1, 50);

   // initial weight (generate random initial weight)
   double sample[50] = { 0.0386, 0.0377, 0.1457, 0.0005, 0.0203, 0.0282, 0.051, 0.2395, 0.1467, 0.1359,
       0.0861, 0.0734, 0.1024, 0.0349, 0.0602, 0.0376, 0.0108, 0.0905, 0.0539, 0.1539,
       0.1113, 0.0153, 0.0675, 0.0605, 0.1976, -0.0118, -0.0346, -0.0102, -0.0104, -0.0309,
       -0.1296, -0.0767, -0.0855, -0.1212, -0.0424, -0.0000, -0.0525, -0.0048, -0.0376, -0.0083,
       -0.0042, -0.0004, -0.0149, -0.0518, -0.0725, -0.0281, -0.1096, -0.0514, -0.004, -0.0066 };

   // generates an initial vector for our invariant distribution 
   double** E0 = Array(1, 50); 
   for (int i = 0; i < 50; ++i) {
       E0[1][i] = sample[i];
   }
   
   // generate our neighbor weight vector
   int i1 = MTUniform(1) * 50; // index of element we will reduce weight
   int i2 = MTUniform(1) * 50; // index of element we will add weight

   E0[1][i1] = E0[1][i1] - 0.0001;
   E0[1][i2] = E0[1][i2] + 0.0001;
   
   // calculate the portoflio variance 
   Multiply(Multiply(E0, V), Transpose(E0));

   // if neighbor state is lower than prior, use initial


   // if neighbor is larger than generate 


   // problem 2.


   // problem 3.


   // problem 4.



   // and ends here. 
   /////////////////////////////////////////////////////////////////////////////

   // Report the best-found portfolio and its variance here.

   // Pause so the execution window does not close.
   Exit ();

}



////////////////////////////////////////////////////////////////////////////////
// Allocate space for and read in covariance matrix and stock tickers.
////////////////////////////////////////////////////////////////////////////////
void GetData () {

   int i, j;
   double V0;
   char input[100];
   FILE *fp;

   // Allocate array space.
   V = Array (50, 50);

   // Allocate space to hold ticker names; "ticker" is a global variable.
   ticker = (char **) calloc (50+1, sizeof (char *));
   for (i = 0; i <= 50; i++) {
      ticker[i] = (char *) calloc (10, sizeof (char));
   }

   // Read in the data.
   fp = fopen ("S&P50Covariances.txt", "r");
   for (i = 1; i <= 50; i++) {

      // Read in stock i's covariance data.

      // Name of the stock ticker.
      fgets (ticker[i], 9, fp);

      // The 50 covariances for stock "i".
      for (j = 1; j <= 50; j++) {

         // Read in V[i][j].
         fgets (input, 99, fp);
         sscanf (input, "%lf", &V0);

         // Put data into the V array.
         V[i][j] = V0;

      }

   }
   fclose (fp);

   return;

}



