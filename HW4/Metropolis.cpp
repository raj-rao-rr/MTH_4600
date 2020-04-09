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

// Calucates the Square of an array 
void SquareArray(double**& arr) {
    double row = arr[0][0], col = arr[0][1];
    int length;

    // checks to see if this is a vertical or horizontal array
    if (row < col) {
        length = col;
    }
    else {
        length = row;
    }

    // modifies the original array in memory 
    for (int i = 1; i <= length; ++i) {
        if (row < col) {
            arr[1][i] = arr[1][i] * arr[1][i];
        }
        else {
            arr[i][1] = arr[i][1] * arr[i][1];
        }
    }
}

// Calculates the Mean Squared Error 
double MSE(double** &actual, double** &pred, double n) {
    // takes the difference between the actual and simualted weights
    double **y1 = Add(Transpose(actual), ScalarMultiple(-1.0, pred));
    // squares the difference before scaling down the terms
    SquareArray(y1);
    double **y2 = ScalarMultiple(1.0 / n, y1);
    double sum = 0.0;

    // computes the sum of reduced terms 
    for (int i = 1; i <= n; ++i) {
        sum += y2[1][i];
    }

    return sum;
}


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
   // variable initialization 
   double **I, **actual_mv, **EX, **c;
   double** e = Array(1, 50);
   double rho, U, var1, var2, error, T = 0.35;

   // define the e-array equal to [1, 1, ... , 1]^T
   for (int i = 1; i <= 50; ++i) {
       e[1][i] = 1;
   }

   // compute the inverted matrix V^{-1}
   I = Invert(V);

   // compute the value of c = e^T*(V^{-1})*e
   c = Multiply(Multiply(e,I), Transpose(e));

   // the actual minimum variance portfolio for the provided covariance matrix 
   actual_mv = Multiply(ScalarMultiple((1/c[1][1]), I), Transpose(e));

   printf("The actual minimum variance portoflio is\n");
   for (int i = 1; i <= 50; ++i) {
       printf("%8.4f\n", actual_mv[i][1]);
   }

   // initial weight (generate random initial weight)
   double sample[50] = { 0.0386, 0.0377, 0.1457, 0.0005, 0.0203, 0.0282, 0.051, 0.2395, 0.1467, 0.1359,
       0.0861, 0.0734, 0.1024, 0.0349, 0.0602, 0.0376, 0.0108, 0.0905, 0.0539, 0.1539,
       0.1113, 0.0153, 0.0675, 0.0605, 0.1976, -0.0118, -0.0346, -0.0102, -0.0104, -0.0309,
       -0.1296, -0.0767, -0.0855, -0.1212, -0.0424, -0.0000, -0.0525, -0.0048, -0.0376, -0.0083,
       -0.0042, -0.0004, -0.0149, -0.0518, -0.0725, -0.0281, -0.1096, -0.0514, -0.004, -0.0066 };

   // generates an initial vector for our invariant distribution 
   double** E0 = Array(1, 50); 
   for (int i = 1; i <= 50; ++i) {
       E0[1][i] = sample[i-1];
   }

   // MONTE CARLO SIMULATION
   for (int j = 1; j <= 10; ++j) {

       // generate our neighbor weight vector
       EX = Copy(E0);
       //double i1 = MTUniform(seed) ; // index of element to reduce weight
       //double i2 = MTUniform(seed) ; // index of element to add weight
       printf("i1 is %8.4f vs i2 is %8.4f\n", MTUniform(seed), MTUniform(seed));
       //// calculate the portoflio variance 
       //var1 = Multiply(Multiply(E0, V), Transpose(E0))[1][1];

       //EX[1][i1] = EX[1][i1] - 0.0001;
       //EX[1][i2] = EX[1][i2] + 0.0001;
       //var2 = Multiply(Multiply(EX, V), Transpose(EX))[1][1];

       //printf("Var1 %8.4f vs Var2 %8.4f\n", var1, var2);

       //// if neighbor state is lower than prior, use the new state as our weight
       //if (var2 < var1) {
       //    E0 = EX;
       //}
       //// if neighbor is larger than prior, use accept/reject scheme 
       //else {
       //    U = MTUniform(seed); // generate a fresh uniform  
       //    rho = exp(-(var2 - var1) / T);

       //    // if U < our rho we modify our weight to the new neighbor 
       //    if (U < rho) {
       //        E0 = EX;
       //    }
       //}

       /*if (j % 100000 == 0) {
           printf("At sim %8.4d with variance %8.4f\n", j, var1);
       }*/
   }

   //// Output the invariant distribution 
   //printf("The calculated minimum variance portoflio is\n");
   //for (int i = 1; i <= 50; ++i) {
   //    printf("%8.4f\n", E0[1][i]);
   //}
   // 
   //// compute the error of the projection
   //error = MSE(actual_mv, E0, 50.0);
   //printf("Our mean squared error is % 8.4f\n", error);

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



