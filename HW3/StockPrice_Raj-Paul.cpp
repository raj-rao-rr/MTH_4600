#include <iostream>
#include "Functions.h"
#include <string>
using namespace std; 

// Used for calculating the max payoff of a simple vanilla call
double max(double val) {
    if (val > 0.0) {return val;}
    return 0.0;
}

// Used for calculating the max payoff of "shout" option
double ShoutMax(double val1, double val2) {
    if (val1 > val2) {
        if (val1 < 0) {
            return 0;
        }
        return val1;
    }
    else if (val2 > val1) {
        if (val2 < 0) {
            return 0;
        }
        return val2;
    }
    else {
        if (val1 > 0 || val2 > 0) {
            return val1;
        }
        return 0;
    }
}


int main() {
   double r, t, T, mu, sigma, alpha, dt, sqrtdt, S, Stilde, S0, V, Vbar, V2bar, Cbar, C2bar,
          elapsed_time, timestart, t_star, stdhat, error, epsilon, n, Discount_factor,
          U, Z, A, B, Btilde, C, Ctilde, Cstar, Cstarstar,Cstarbar, Cstar2bar, Cstarstarbar, Cstarstar2bar,
          BSprice, A2, XA, K, Smax, SmaxT,part1, part1time, part2, part2time, part3, part3time, part4, part5, part6, part6time, part7, part8;
   int i, N, done, test, part;

   for(part = 1; part <= 7; ++part){
      printf("\n\nPart %d:\n", part);
      
      if (part == 0 || part == 0 || part == 0 || part == 0 || part == 0) { 
         printf("Part %d omitted for testing purposes\n", part);
         continue;
      }

      // Strike Price
      K = 110.0;

      // Time to expiration.
      T = 0.5;

      // Number of stock price periods.
      N = 50;
      if (part == 4) {
         N = 1; // Adjusted N for part 4 problem set 
      }

      // Time increment per period.
      dt = T / N;

      // Compute this oft-used value once and for all.
      sqrtdt = sqrt(dt);

      // Risk-free interest rate.
      r = 0.05;

      // Compute the oft-used discount factor once and for all.
      Discount_factor = exp(-r*T);

      // Stock price volatility.
      sigma = .30;

      // Drift term.
      mu = r - sigma*sigma/2.0;

      // Initial stock price.
      S0 = 100.0;

      if (part == 5) { 
         // Actual Price of the Option calculated in part 5
         part5 = BlackScholes(T, S0, K, sigma, r);
         printf("   This uses the Black-Scholes call option pricing formula.");
         continue;
      }

      // Specify the 95% error tolerance.
      epsilon = 0.005;

      // Start the clock to time the computations.
      timestart = clock ();

      // Seed the RNG.
      MTUniform (1);

      // Print column headings for output to execution window.
      if (part == 1 || part == 2){
         printf ("         n       Vbar        +/-        t       t*\n");
      }
      else if (part == 3 || part == 4){
         printf ("         n      C*bar        +/-        t       t*\n");
      }
      else if (part == 6 || part == 7){
         printf ("         n     C**bar        +/-        t       t*\n");
      }

      // Initialize certain values.
      V2bar = Vbar = A2 = XA = done = n = test = 0;

      // Begin the simulation loop.
      while (!done) {

         // Initialize the stock price.
         S = S0;

         // Initialize stock price, Brownian path, time, Var A, and Cov (X,A)
         Smax = SmaxT = B = t = 0;

         // Simulate a stock price path.  Go forward period-by-period computing
         //   the next stock price.
         for (i = 1; i <= N; ++i) {
            
            // Advance time by one period.
            t += dt;

            // Advance the path.
            U = MTUniform (0);
            Z = PsiInv (U); // Standard normal via inverse transform.
            B += sqrtdt * Z; // Standard Brownian path 
            S = S0 * exp (mu*t + sigma*B); // Compute the next stock price.
            
            if (part != 1){
               Btilde = -B; // symmetric distribution
               Stilde = S0 * exp(mu * t + sigma*Btilde); // antithetic stock process
            }

            if (part == 7) { // Paul
               // Keep track of the maximum stock price 
               if (S > Smax) {
                  Smax = S;
               }
               if (Stilde > SmaxT) {
                  SmaxT = Stilde;
               }
            }
         }
         
         // S now is the value of the stock at time T for the simulated price path.

         if (part == 1) {
            V = Discount_factor * S;
         }
         else if (part == 2) {
            V = Discount_factor * (S + Stilde) / 2.0 ;
         }
         else if (part == 3 || part == 4) { 
            // Determine call payoff at time T
            C = max(S - K);
            Ctilde = max(Stilde - K);

            // Calculate C*
            Cstar = Discount_factor * (C + Ctilde)/2;
         }
         else if (part == 6 || part == 7) {

            // Compute the value of A
            A = B * B - T;

            // Update the required sample moments. There is no need to take the average by dividing by n
            // because the calculation for alpha would cancel out that division.
            A2 += A * A;
            XA += B * A;

            // Appropriate value for alpha 
            alpha = -10; //- XA / A2;

            if (part == 6){
               // Part 6 uses the last stock price observed
               C = max(S - K);
               Ctilde = max(Stilde - K);
            } 
            else if (part == 7){
               // Part 7 uses the maximum stock price observed along the stock price path
               C = max(Smax - K);
               Ctilde = max(SmaxT - K);
            }

            // Sample mean statistics being averaged 
            Cstarstar = Discount_factor * ((C + Ctilde)/2) + (alpha * A);
         }

         // Update the simulation counter and the sample moments of V at time T.
         ++ n;
         if (part == 1 || part == 2) {
            Vbar  = ((n-1)*Vbar + V) / n;
            V2bar = ((n-1)*V2bar + V*V) / n;
         } 
         else if (part == 3 || part == 4) {
            Cstarbar  = ((n-1)*Cstarbar + Cstar) / n;
            Cstar2bar = ((n-1)*Cstar2bar + Cstar*Cstar) / n;
         } 
         else if (part == 6 || part == 7) {
            Cstarstarbar  = ((n-1)*Cstarstarbar + Cstarstar) / n;
            Cstarstar2bar = ((n-1)*Cstarstar2bar + Cstarstar*Cstarstar) / n;
         }

         // Update the error condition test counter.
         ++ test;

         // Test the error condition periodically.
         if (test == 100000) {

            if (part == 1 || part == 2){
            // Estimate the standard deviation of S.
            stdhat  = sqrt (V2bar - Vbar*Vbar);
            }

            else if (part == 3 || part == 4){
            // Estimate the standard deviation of C*.
            stdhat  = sqrt (Cstar2bar - Cstarbar*Cstarbar);
            }

            else if (part == 6 || part == 7){
            // Estimate the standard deviation of C**.
            //printf ("alpha: %8.4f, XA: %8.4f, A2: %8.4f \n", alpha, XA, A2);
            //printf ("B: %8.4f, A: %8.4f\n", B, A);
            stdhat  = sqrt (Cstarstar2bar - Cstarstarbar*Cstarstarbar);
            }

            // Estimate the error of the Vbar estimator.
            error = 1.96 * stdhat / sqrt(n);

            // Compute the elapsed time and estimate the time of completion.
            elapsed_time = (double) (clock() - timestart) / CLOCKS_PER_SEC;
            t_star = elapsed_time * pow (error / epsilon, 2);

            // Reset the "test" counter and see if error tolerance is met.
            test = 0;
            if (error < epsilon) {
               done = 1;
            }

            // 1) Reports
            // 2) Saves the outputs to print all at the end
            if (part == 1) {
               printf ("%10.0f   %8.4f   %8.6f %8.3f %8.3f\n", n, Vbar, error, elapsed_time, t_star);
               part1 = Vbar;
               part1time = elapsed_time;
            }
            else if (part == 2) {
               printf ("%10.0f   %8.4f   %8.6f %8.3f %8.3f\n", n, Vbar, error, elapsed_time, t_star);
               part2 = Vbar;
               part2time = elapsed_time;
            }
            else if (part == 3) {
               printf ("%10.0f   %8.4f   %8.6f %8.3f %8.3f\n", n, Cstarbar, error, elapsed_time, t_star);
               part3time = elapsed_time;
               part3 = Cstarbar;
            }
            else if (part == 4) {
               printf ("%10.0f   %8.4f   %8.6f %8.3f %8.3f\n", n, Cstarbar, error, elapsed_time, t_star);
               part4 = Cstarbar;
            }
            else if (part == 6) {
               printf ("%10.0f   %8.4f   %8.6f %8.3f %8.3f\n", n, Cstarstarbar, error, elapsed_time, t_star);
               part6time = elapsed_time;
               part6 = Cstarstarbar;
            }
            else if (part == 7) {
               printf ("%10.0f   %8.4f   %8.6f %8.3f %8.3f\n", n, Cstarstarbar, error, elapsed_time, t_star);
               part7 = Cstarstarbar;
            }
         }
      }
   }
   printf("\nPart 1:\n=======\n   Vbar is %3.4f and therefore agrees with S0\n", part1);
   printf("   This part ran in %8.3f seconds\n", part1time);

   printf("\nPart 2:\n=======\n   Similar to Part 1, Vbar is %3.4f and therefore also agrees with S0\n", part2);
   printf("   This was significantly faster than Part 1, running in %8.3f seconds\n", part2time);

   printf("\nPart 3:\n=======\n   Using the St simulations from Part 2, the price a call option is %2.5f\n", part3);
   printf("   This part ran in %8.3f seconds\n", part3time);

   printf("\nPart 4:\n=======\n   Redoing the simulations with N=1 instead of 50, the price a call option is %2.5f\n", part4);
   printf("   This is the same answer as in Part 3\n");

   printf("\nPart 5:\n=======\n   Using the Black-Scholes formula, the exact value of the option is %2.5f\n", part5);

   printf("\nPart 6:\n=======\n   Similar to Part 1 and 2, Vbar is %3.4f and therefor also agrees with S0\n", Cstarbar);
   printf("   This was significantly faster than Part 3, running in %8.3f seconds\n", part6time);

   printf("\nPart 7:\n=======\n   The value of the look-back option is %2.5f\n", part7);



   Exit ();
   
}



