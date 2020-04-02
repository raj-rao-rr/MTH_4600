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



