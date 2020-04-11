////////////////////////////////////////////////////////////////////////////////
// This is starter code for Project 5.
////////////////////////////////////////////////////////////////////////////////


#include "Functions.h"

// This function is found below.
void GetData ();

// "ticker" is a global variable.
char **ticker;

// The covariance matrix is a global variable.
double **V;

// Computes the Variance
double Variance(double** &arr, double** &cov) {
    // computes the expression wCw^T, where w - weights and C - covariance 
    double** b = Transpose(arr);
    double** a = Multiply(b, cov);
    //double** b = Transpose(arr);
    return Multiply(a, arr)[1][1];
}


int main () {

   int i, k, l, seed, n, n_min, Stock_A, Stock_B, accept;
   double **Weights, **BestWeights, **e, **I, **actual_mv;
   double Var, Var_min, DeltaVar, Var_PreChange, Var_PostChange, U, T, c;

   // Read in and display the tickers and covariance matrix.
   GetData ();

   // Show the 50 tickers.
   for (i = 1; i <= 50; i++) {
      printf (ticker[i]);
   }
   Pause ();

   // Show the covariance matrix.
   Show (V);

   // define the e-array equal to [1, 1, ... , 1]^T
   e = Array(1, 50);
   for (int i = 1; i <= 50; ++i) {
      e[1][i] = 1;
   }

   // compute the inverted matrix V^{-1}
   I = Invert(V);

   // compute the value of c = e^T*(V^{-1})*e
   c = Multiply(Multiply(e,I), Transpose(e))[1][1];

   // the actual minimum variance portfolio for the provided covariance matrix 
   actual_mv = Multiply(ScalarMultiple((1/c), I), Transpose(e));

   printf("The actual minimum variance portoflio is\n");
   for (int i = 1; i <= 50; ++i) {
      printf("%8.4f\n", actual_mv[i][1]);
   }

   // Seed the RNG.
   seed = 0; // GetInteger ("RNG Seed?... ");
   MTUniform (seed);

   // Temperature.
   T = 0.000000001; //GetDouble ("Temperature? (Good choice: 0.000000005)... ");

   /////////////////////////////////////////////////////////////////////////////
   // Your Metropolis algorithm starts here...

   // Intitialize Weight Array, with all 50 stocks having equal weighting (2% = 200 cents)
   Weights = Array(50,1);
   for (i = 1; i <= 50; i++){
      Weights[i][1] = 0.02;
   }

   Var = 0;
   Var = Variance(Weights, V);

   // Initialize the minimal variance and where it occurs.
   Var_min = Var;
   n_min = 0;

   // Record the initial route as the best-so-far.
   BestWeights = Array(50,1);
   for (i = 1; i <= 50; i++) {
      BestWeights[i, 1] = Weights[i, 1];
   }

   // Simulate the Markov chain 100 million steps.
   for (n = 1; n <= 10000000; n++) {

      // Determine the proposed transition...
      // Randomly choose a "neighbor" of the current route; use accept/reject.
      // First pick a pair of stocks (Stock_A, Stock_B), uniformly from {1,...,50} x {1,...,50} 
      // with Stock_A != Stock_B
      while (1) {

         // Pick Stock_A_Weight and Stock_B_Weight independently and uniformly from {1,...,50}.
         U = MTUniform (0);
         Stock_A = 1 + int ((50) * U);
         U = MTUniform (0);
         Stock_B = 1 + int ((50) * U);

         // See if they are acceptable, i.e, if they satisfy (1) and (2) above.
         if (Stock_A != Stock_B) {
            break;
         }

      }

      Var_PreChange = Var;

      // Decrease Stock_A's weight by 0.0001 (1 cent) and 
      // increase Stock_B's weight by 0.0001 (1 cent) 
      Weights[Stock_A][1] -= 0.0001;
      Weights[Stock_B][1] += 0.0001;

      Var_PostChange = Variance(Weights, V);

      // Compute the change in variance associated with reversing that weight change
      DeltaVar =  Var_PostChange - Var_min;

      // See if the proposed transition is accepted.
      accept = 0;

      if (DeltaVar <= 0) {
         accept = 1;
      }

      else if (T > 0) {
         if (MTUniform (0) <= exp (-DeltaVar / T)) {
            accept = 1;
         }
      }

      // K:
      // Effect the transition. If the new energy is a best-so-far, record the
      //    route (in best[*]) and the new minimal energy (E_min).
      if (accept) {

         // Reverse the part of route "c" from i to j.
         Var = Var_PostChange;

         // Record data for the best-so-far, if appropriate.
         if (Var_PostChange < Var_min) {
            Var_min = Var_PostChange;
            n_min = n;
            for (l = 1; l <= 50; l++) {
               BestWeights[l][1] = Weights[l][1];
            }
         }
      } // End of "if" statement.
      else {
         Weights[Stock_A][1] += 0.00001;
         Weights[Stock_B][1] -= 0.00001;
         Var = Var_PreChange;
      }

      if (n % 100000 == 0) {
         printf("At sim %8.4d -> variance is %8.10f\n", n, Var);
         printf("                -> Var_min  is %8.10f\n", Var_min);
      }

   }

   printf("Portoflio weights is\n");
   printf("Metropolis - LAST        Metropolis - BEST        Actual       Difference\n");
   for (int i = 1; i <= 50; ++i) {
      printf("%8.6f                 %8.6f                 %8.6f        %8.6f\n", Weights[i][1], BestWeights[i][1], actual_mv[i][1], BestWeights[i][1] - actual_mv[i][1]);
   }

   printf("\n\nCent allocations\n");
   printf("Metropolis - LAST        Metropolis - BEST        Actual       Difference\n");
   for (int i = 1; i <= 50; ++i) {
      printf("%8.0f                 %8.0f              %8.0f       %8.0f\n", Weights[i][1]*100000, BestWeights[i][1]*100000, actual_mv[i][1]*100000, BestWeights[i][1]*100000 - actual_mv[i][1]*100000);
   }

   printf("\n\nPorfolio Vars\n");
   printf("Metropolis - LAST        Metropolis - BEST        Actual       Difference\n");
   printf("%8.10f                %8.10f             %8.10f        %8.10f\n", Variance(Weights, V), Variance(BestWeights, V), Variance(actual_mv, V), Variance(BestWeights, V) - Variance(actual_mv, V));
   
   
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
   fp = fopen ("/Users/Paul/Desktop/MTH 4600/Visual Studio Code Folder/HW4/S&P50Covariances.txt", "r");
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



