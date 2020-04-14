////////////////////////////////////////////////////////////////////////////////
// This is starter code for Project 5.
////////////////////////////////////////////////////////////////////////////////

// This function is found below.
void GetData ();
void SquareArray(double**&);
double Variance(double** , double**);
double MSE(double**&, double**&, double);
int sFlag(double**&, int);

// "ticker" is a global variable.
char **ticker;

// The covariance matrix is a global variable.
double **V;

#include "Functions.h"
#include <iostream>


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

   //////////////////////////////////////////////////////////////////
   // problem 1. 
   //////////////////////////////////////////////////////////////////

   //// variable initialization 
   //double **I, **actual_mv;
   //double **e = Array(1, 50);
   //double rho, c, U, var1, var2, error, deltaVar, T = 0.000001;
   //int i1, i2, seed=0;

   //// define the e-array equal to [1, 1, ... , 1]^T
   //for (int i = 1; i <= 50; ++i) {
   //    e[1][i] = 1;
   //}

   //// compute the inverted matrix V^{-1}
   //I = Invert(V);

   //// compute the value of c = e^T*(V^{-1})*e
   //c = Multiply(Multiply(e,I), Transpose(e))[1][1];

   //// the actual minimum variance portfolio for the provided covariance matrix 
   //actual_mv = Multiply(ScalarMultiple((1/c), I), Transpose(e));

   //printf("The actual minimum variance portoflio is\n");
   //for (int i = 1; i <= 50; ++i) {
   //    printf("%8.4f\n", actual_mv[i][1]);
   //}
   //printf("The minimum variance reached is %8.8f\n", Variance(actual_mv, V));
   //printf("-------------------------------------------------------------\n");

   //// generates an initial vector for our invariant distribution 
   //double** E0 = Array(50, 1); 
   //for (int i = 1; i <= 50; ++i) {
   //    E0[i][1] = 0.0200;
   //}

   //double** EX = Array(50, 1);
   //for (int i = 1; i <= 50; ++i) {
   //    EX[i][1] = E0[i][1];
   //}

   //// calculate the portfolio variance for original weight
   //var1 = Variance(E0, V);

   //printf("Simulation begins...\n");
   //// MONTE CARLO SIMULATION
   //for (int j = 1; j <= 1100000; ++j) {

   //    // selects our stock weights to swap 
   //    while (1) {
   //        // Pick two numbers independently and uniformly from {1,...,50}.
   //        i1 = 1 + int(MTUniform(seed) * 50); // index of element to reduce weight
   //        i2 = 1 + int(MTUniform(seed) * 50); // index of element to add weight
  
   //        // See if they are acceptable, i.e, if they satisfy (1) and (2) above.
   //        if (i1 != i2) {
   //            break;
   //        }
   //    }
   //   
   //    // building the neighbor to the state
   //    EX[i1][1] -= 0.0001;
   //    EX[i2][1] += 0.0001;

   //    // calculate the portfolio variance for weight 2
   //    var2 = Variance(EX, V);

   //    // compute the change in variances
   //    deltaVar = var2 - var1;

   //    // if neighbor state is lower than prior, use the new state as our weight
   //    if (deltaVar <= 0) {
   //        var1 += deltaVar;
   //    }
   //    // if neighbor is larger than prior, use accept/reject scheme 
   //    else {
   //        U = MTUniform(seed); // generate a fresh uniform  
   //        rho = exp(-deltaVar / T);

   //        // if U < our rho we modify our weight to the new neighbor 
   //        if (U <= rho) {
   //            var1 += deltaVar;
   //        }
   //        else {
   //            // returning back to the previous state 
   //            EX[i1][1] += 0.0001;
   //            EX[i2][1] -= 0.0001;
   //        }
   //    }

   //    if (j % 100000 == 0) {
   //        printf("At sim %8.4d -> variance is %8.8f\n", j, var1);
   //    }
   //}

   //// Report the best-found portfolio and its variance here.
   //printf("The calculated minimum variance portoflio is\n");
   //for (int i = 1; i <= 50; ++i) {
   //    printf("%8.4f\n", EX[i][1]);
   //}
   //printf("The minimum variance reached is %8.8f\n", var1);
   // 
   //// compute the error of the projection
   //error = MSE(actual_mv, EX, 50.0);
   //printf("Our mean squared error is % 8.8f\n", error);
   


   //////////////////////////////////////////////////////////////////
   // problem 2.
   //////////////////////////////////////////////////////////////////

   //// variable initialization 
   //double rho, U, var1, var2, error, deltaVar, T = 0.000001;
   //int i1, i2, flag, seed=0;

   //// generates an initial vector for our invariant distribution 
   //double** E0 = Array(50, 1); 
   //for (int i = 1; i <= 50; ++i) {
   //    E0[i][1] = 0.0200;
   //}

   //double** EX = Array(50, 1);
   //for (int i = 1; i <= 50; ++i) {
   //    EX[i][1] = E0[i][1];
   //}

   //// calculate the portfolio variance for original weight
   //var1 = Variance(E0, V);

   //printf("Simulation begins...\n");
   //// MONTE CARLO SIMULATION
   //for (int j = 1; j <= 2000000; ++j) {

   //    // selects our stock weights to swap 
   //    while (1) {
   //        // Pick two numbers independently and uniformly from {1,...,50}.
   //        i1 = 1 + int(MTUniform(seed) * 50); // index of element to reduce weight
   //        i2 = 1 + int(MTUniform(seed) * 50); // index of element to add weight
  
   //        // See if they are acceptable, i.e, if they satisfy (1) and (2) above.
   //        if (i1 != i2) {
   //            break;
   //        }
   //    }
   //   
   //    // building the neighbor to the state
   //    EX[i1][1] -= 0.0001;
   //    EX[i2][1] += 0.0001;

   //    // Check to see if there is a short position in the weights
   //    flag = sFlag(EX, 50);

   //    if (flag == 0) {
   //        // calculate the portfolio variance for weight 2
   //        var2 = Variance(EX, V);
   //    }
   //    else {
   //        var2 = 1000;
   //    }

   //    // compute the change in variances
   //    deltaVar = var2 - var1;

   //    // if neighbor state is lower than prior, use the new state as our weight
   //    if (deltaVar <= 0) {
   //        var1 += deltaVar;
   //    }
   //    // if neighbor is larger than prior, use accept/reject scheme 
   //    else {
   //        U = MTUniform(seed); // generate a fresh uniform  
   //        rho = exp(-deltaVar / T);

   //        // if U < our rho we modify our weight to the new neighbor 
   //        if (U <= rho) {
   //            var1 += deltaVar;
   //        }
   //        else {
   //            // returning back to the previous state 
   //            EX[i1][1] += 0.0001;
   //            EX[i2][1] -= 0.0001;
   //        }
   //    }

   //    if (j % 100000 == 0) {
   //        printf("At sim %8.4d -> variance is %8.8f\n", j, var1);
   //    }
   //}

   //// Report the best-found portfolio and its variance here.
   //printf("The calculated minimum variance portoflio is\n");
   //for (int i = 1; i <= 50; ++i) {
   //    printf("%8.4f\n", EX[i][1]);
   //}
   //printf("The minimum variance reached is %8.8f\n", var1);



   //////////////////////////////////////////////////////////////////
   // problem 3.
   //////////////////////////////////////////////////////////////////

   //// variable initialization 
   //double rho, U, var1, var2, error, deltaVar, wt = 1.0, mod = 0.0, num = 1.0, T = 0.02;
   //int i1, seed = 0;

   //// generates an initial vector for our invariant distribution 
   //double** E0 = Array(50, 1);
   //E0[50][1] = wt;

   //double** EX = Array(50, 1);
   //for (int i = 1; i <= 50; ++i) {
   //    EX[i][1] = E0[i][1];
   //}

   //// calculate the portfolio variance for original weight
   //var1 = Variance(E0, V);

   //printf("Simulation begins...\n");
   //// MONTE CARLO SIMULATION
   //for (int j = 1; j <= 1100000; ++j) {

   //    // Pick two numbers independently and uniformly from {1,...,50}.
   //    i1 = 1 + int(MTUniform(seed) * 50); // index of element to reduce weight

   //    // performs a "flip" in defining the neighbors for each state
   //    if (EX[i1][1] > 0) {
   //        EX[i1][1] = 0.0;
   //        mod = -1;
   //        num += mod;
   //    }
   //    else {
   //        EX[i1][1] = 1.0;
   //        mod = 1;
   //        num += mod;
   //    }

   //    // if their are no stocks present we assign no wt to the portfolio
   //    if (num != 0) {
   //        // simple portfolio even distribution
   //        wt = 1.0 / num;
   //    }
   //    else {
   //        wt = 0.0;
   //    }

   //    // reassign the weights of the vector evenly
   //    for (int i = 1; i <= 50; ++i) {
   //        // if there is a postive weight we reasign it 
   //        if (EX[i][1] > 0) {
   //            EX[i][1] = wt;
   //        }
   //    }

   //    if (num != 0) {
   //        // calculate the portfolio variance for weight 2
   //        var2 = Variance(EX, V);
   //    }
   //    else {
   //        var2 = 1000;
   //    }

   //    // compute the change in variances
   //    deltaVar = var2 - var1;

   //    // if neighbor state is lower than prior, use the new state as our weight
   //    if (deltaVar <= 0) {
   //        var1 += deltaVar;
   //    }
   //    // if neighbor is larger than prior, use accept/reject scheme 
   //    else {
   //        U = MTUniform(seed); // generate a fresh uniform  
   //        rho = exp(-deltaVar / T);

   //        // if U < our rho we modify our weight to the new neighbor 
   //        if (U <= rho) {
   //            var1 += deltaVar;
   //        }
   //        else {
   //            // returning back to the previous state (number of stocks and weights)
   //            num -= mod;
   //            wt = 1.0 / num;

   //            // if the modification was to add a stock, we reduce the weight to 0.0
   //            if (mod > 0) {
   //                EX[i1][1] = 0.0;
   //            }
   //            // if the modification was to remove a stock, we add the weight back 
   //            else {
   //                EX[i1][1] = 1.0;
   //            }
   //            
   //            // reassign the weights of the vector evenly
   //            for (int i = 1; i <= 50; ++i) {
   //                // if there is a postive weight we reasign it 
   //                if (EX[i][1] > 0) {
   //                    EX[i][1] = wt;
   //                }
   //            }

   //        }
   //    }

   //    if (j % 100000 == 0) {
   //        printf("At sim %8.4d -> variance is %8.8f\n", j, var1);
   //    }
   //}

   //// Report the best-found portfolio and its variance here.
   //printf("The calculated minimum variance portoflio is\n");
   //for (int i = 1; i <= 50; ++i) {
   //    printf("%8.4f\n", EX[i][1]);
   //}
   //printf("The minimum variance reached is %8.8f\n", var1);


   //////////////////////////////////////////////////////////////////
   // problem 4.
   //////////////////////////////////////////////////////////////////

   // variable initialization 
   double rho, U, var1, var2, var3, error, deltaVar, wt = 1.0, mod = 0.0, num = 1.0, T = 0.00;
   int i1, seed = 0;

   // generates an initial vector for our invariant distribution 
   double** E0 = Array(50, 1);
   E0[50][1] = wt;

   double** EX = Array(50, 1);
   for (int t = 1; t <= 50; ++t) {
       EX[t][1] = E0[t][1];
   }

   // calculate the portfolio variance for original weight
   var1 = Variance(E0, V);

   printf("Simulation begins...\n");
   // MONTE CARLO SIMULATION
   for (int j = 1; j <= 500000; ++j) {

       // Pick two numbers independently and uniformly from {1,...,50}.
       i1 = 1 + int(MTUniform(seed) * 50); // index of element to reduce weight

       // performs a "flip" in defining the neighbors for each state
       if (EX[i1][1] > 0) {
           EX[i1][1] = 0.0;
           mod = -1;
           num += mod;
       }
       else {
           EX[i1][1] = 1.0;
           mod = 1;
           num += mod;
       }

       // if their are no stocks present we assign no wt to the portfolio
       if (num != 0) {
           // simple portfolio even distribution
           wt = 1.0 / num;
       }
       else {
           wt = 0.0;
       }

       // reassign the weights of the vector evenly
       for (int c = 1; c <= 50; ++c) {
           // if there is a postive weight we reasign it 
           if (EX[c][1] > 0) {
               EX[c][1] = wt;
           }
       }

       if (num != 0) {
           // calculate the portfolio variance for weight 2
           var2 = Variance(EX, V);
       }
       else {
           var2 = 1000;
       }

       // compute the change in variances
       deltaVar = var2 - var1;

       // if neighbor state is lower than prior, use the new state as our weight
       if (deltaVar <= 0) {
           var1 += deltaVar;
       }
       // if neighbor is larger than prior, use accept/reject scheme 
       else {
           U = MTUniform(seed); // generate a fresh uniform  
           rho = exp(-deltaVar / T);

           // if U < our rho we modify our weight to the new neighbor 
           if (U <= rho) {
               var1 += deltaVar;
           }
           else {
               // returning back to the previous state (number of stocks and weights)
               num -= mod;
               wt = 1.0 / num;

               // if the modification was to add a stock, we reduce the weight to 0.0
               if (mod > 0) {
                   EX[i1][1] = 0.0;
               }
               // if the modification was to remove a stock, we add the weight back 
               else {
                   EX[i1][1] = 1.0;
               }
               
               // reassign the weights of the vector evenly
               for (int v = 1; v <= 50; ++v) {
                   // if there is a postive weight we reasign it 
                   if (EX[v][1] > 0) {
                       EX[v][1] = wt;
                   }
               }

           }
       }

       if (j % 100000 == 0) {
           printf("At sim %8.4d -> variance is %8.8f...\n", j, var1);
       }
   }

   delete E0;
   // Verify that the state constructed is in fact a stable state 
   double** EZ = Array(50, 1);
   for (int w = 1; w <= 50; ++w) {
       EZ[w][1] = EX[w][1];
   }

   printf("Checking the stability of the state...\n");
   // count the amount of times the var is best
   int count = 0;

   // itterate throuhg each neighboring state
   for (int j = 1; j <= 50; ++j) {

       // performs a "flip" in defining the neighbors for each state
       if (EZ[j][1] > 0) {
           EZ[j][1] = 0.0;
           mod = -1;
           num += mod;
       }
       else {
           EZ[j][1] = 1.0;
           mod = 1;
           num += mod;
       }

       // if their are no stocks present we assign no wt to the portfolio
       if (num != 0) {
           // simple portfolio even distribution
           wt = 1.0 / num;
       }
       else {
           wt = 0.0;
       }

       // reassign the weights of the vector evenly
       for (int i = 1; i <= 50; ++i) {
           // if there is a postive weight we reasign it 
           if (EZ[i][1] > 0) {
               EZ[i][1] = wt;
           }
       }

       if (num != 0) {
           // calculate the portfolio variance for weight 2
           var3 = Variance(EZ, V);
       }
       else {
           var3 = 1000;
       }

       // check to see if the state is stable, checks against neighors
       if (var3 > var1) {
           // counting the number of neighbors that have less variance
           ++count;
           
           // returning back to the previous state (number of stocks and weights)
           num -= mod;
           wt = 1.0 / num;

           // if the modification was to add a stock, we reduce the weight to 0.0
           if (mod > 0) {
               EZ[j][1] = 0.0;
           }
           // if the modification was to remove a stock, we add the weight back 
           else {
               EZ[j][1] = 1.0;
           }

           // reassign the weights of the vector evenly
           for (int i = 1; i <= 50; ++i) {
               // if there is a postive weight we reasign it 
               if (EZ[i][1] > 0) {
                   EZ[i][1] = wt;
               }
           }
       }
       else {
           printf("A neighboring state has violated the condition, variance was %8.6f\n", var3);
           break;
       }

   }

   // if the state is lower than all neighbors it is stable 
   if (count == 50) {
       printf("We have found a stable state, the portfolio is\n");
       printf("Stock  |  Weighting \n");
       printf("-----------------\n");
       // Report the best-found portfolio and its variance here.
       for (int i = 1; i <= 50; ++i) {
           printf("%s %8.4f\n", ticker[i], EX[i][1]);
       }
       printf("The minimum variance reached is %8.8f\n", var1);
   }

   // Pause so the execution window does not close.
   Exit ();

}

////////////////////////////////////////////////////////////////////////////////
// Allocate space for helper functions
////////////////////////////////////////////////////////////////////////////////


// Calucates the Square of an array 
void SquareArray(double**& arr) {
    // modifies the original array in memory 
    for (int i = 1; i <= 50; ++i) {
        arr[i][1] = arr[i][1] * arr[i][1];
    }
}

// Identifies the presence of a Short position, flags with 1 for True
int sFlag(double**& arr, int size) {
    int flag = 0;
    for (int i = 1; i <= size; ++i) {
        if (arr[i][1] < 0.0) {
            flag = 1;
            break;
        }
    }
    return flag;
}

// Computes the Variance of a given set of weights
double Variance(double** arr, double** cov) {
    // computes the expression wCw^T, where w - weights and C - covariance 
    double** b = Transpose(arr);
    double** a = Multiply(b, cov);
    return Multiply(a, arr)[1][1];
}

// Calculates the Mean Squared Error 
double MSE(double**& actual, double**& pred, double n) {
    // takes the difference between the actual and simualted weights
    double** y1 = Add(actual, ScalarMultiple(-1.0, pred));

    // squares the difference before scaling down the terms
    SquareArray(y1);
    double** y2 = ScalarMultiple(1.0 / n, y1);

    double sum = 0.0;
    // computes the sum of reduced terms 
    for (int i = 1; i <= n; ++i) {
        sum += y2[i][1];
    }

    return sum;
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



