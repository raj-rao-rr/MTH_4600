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
#include <iostream>


// Calucates the Square of an array 
void SquareArray(double**& arr) {
    // modifies the original array in memory 
    for (int i = 1; i <= 50; ++i) {
        arr[i][1] = arr[i][1] * arr[i][1];
    }
}

// Identifies the Presence of a Short position, flags with 1 for True
int sFlag(double** arr, int size) {
    int flag = 0;
    for (int i = 1; i < size; ++i) {
        if (arr[i][1] < 0) {
            flag = 1;
            break;
        }
    }
    return flag;
}

// Computes the Variance
double Variance(double** &arr, double** &cov) {
    // computes the expression wCw^T, where w - weights and C - covariance 
    double** b = Transpose(arr);
    double** a = Multiply(b, cov);
    return Multiply(a, arr)[1][1];
}

// Calculates the Mean Squared Error 
double MSE(double** &actual, double** &pred, double n) {
    // takes the difference between the actual and simualted weights
    double **y1 = Add(actual, ScalarMultiple(-1.0, pred));

    // squares the difference before scaling down the terms
    SquareArray(y1);
    double **y2 = ScalarMultiple(1.0 / n, y1);

    double sum = 0.0;
    // computes the sum of reduced terms 
    for (int i = 1; i <= n; ++i) {
        sum += y2[i][1];
    }

    return sum;
}


int main () {

   int i;

   // Read in and display the tickers and covariance matrix.
   GetData ();

   // Show the 50 tickers.
   for (i = 1; i <= 50; ++i) {
      printf (ticker[i]);
   }
   Pause ();

   // Show the covariance matrix.
   Show (V);

   /////////////////////////////////////////////////////////////////////////////
   // Your Metropolis algorithm starts here...

   // problem 1. 
   // variable initialization 
   double **I, **actual_mv;
   double **e = Array(1, 50);
   double rho, c, U, var1, var2, error, T = 0.35;
   int i1, i2;

   // define the e-array equal to [1, 1, ... , 1]^T
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
   printf("The minimum variance reached is %8.4f\n", Variance(actual_mv, V));

   // generates an initial vector for our invariant distribution 
   double **EX,** E0 = Array(50, 1); 
   for (int i = 1; i <= 50; ++i) {
       E0[i][1] = 0.02;
   }

   // calculate the portfolio variance for original weight
   var1 = Variance(E0, V);

   printf("Simulation begins...\n");
   // MONTE CARLO SIMULATION
   for (int j = 1; j <= 100; ++j) {
       std::cout << "Memory location of E0 is \n" << (long) (&E0) << "\n";
       std::cout << "Memory location of EX is \n" << (long) (&EX) << "\n";
       // generate our neighbor weight vector
       EX = E0;

       // selects our stock weights to swap 
       while (1) {
           // Pick two numbers independently and uniformly from {1,...,50}.
           i1 = 1 + (int) (MTUniform(0) * 50); // index of element to reduce weight
           i2 = 1 + (int) (MTUniform(0) * 50); // index of element to add weight

           // See if they are acceptable, i.e, if they satisfy (1) and (2) above.
           if (i1 != i2) {
               break;
           }
       }
      
       // building the neighbor to the state
       EX[i1][1] -= 0.0001;
       EX[i2][1] += 0.0001;

       // calculate the portfolio variance for weight 2
       var2 = Variance(EX, V);

       // if neighbor state is lower than prior, use the new state as our weight
       if (var2 < var1) {
           E0 = EX;
           var1 = var2;
       }
       // if neighbor is larger than prior, use accept/reject scheme 
       else {
           U = MTUniform(0); // generate a fresh uniform  
           rho = exp(-(var2 - var1) / T);

           // if U < our rho we modify our weight to the new neighbor 
           if (U < rho) {
               E0 = EX;
               var1 = var2;
           }
       }

       /*if (j % 10 == 0) {
           printf("At sim %8.4d -> variance is %8.4f\n", j, var1);
       }*/
   }

   // Output the invariant distribution 
   printf("The calculated minimum variance portoflio is\n");
   for (int i = 1; i <= 50; ++i) {
       printf("%8.4f\n", E0[i][1]);
   }
   printf("The minimum variance reached is %8.4f\n", Variance(E0, V));
    
   // compute the error of the projection
   error = MSE(actual_mv, E0, 50.0);
   printf("Our mean squared error is % 8.4f\n", error);

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
   for (i = 1; i <= 50; ++i) {

      // Read in stock i's covariance data.

      // Name of the stock ticker.
      fgets (ticker[i], 9, fp);

      // The 50 covariances for stock "i".
      for (j = 1; j <= 50; ++j) {

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



