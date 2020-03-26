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
    if (val1 > val2) { return val1; }
    else if (val2 > val1) { return val2; }
    else { return 0.0; }
}

int main() {
   double r, t, T, mu, sigma, alpha, alphafound, dt, sqrtdt, S, Stilde, S0, V, Vbar, V2bar, Cbar, C2bar,
          elapsed_time, timestart, t_star, stdhat, error, epsilon, n, Discount_factor,
          U, Z, A, B, Btilde, C, Ctilde, Cstar, Cstarstar,Cstarbar, Cstar2bar, Cstarstarbar, Cstarstar2bar,
          BSprice, A2, XA, K, Smax, SmaxT,part1, part1time, part2, part2time, part3, part3time, part4, part5,
          part6, part6time, part6alpha, part7, part7alpha, part8, X, A2bar, Xbar, XAbar, X2bar, CorA_Cstar;
   int i, N, done, test, part;

   alphafound = 0;

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
         N = 1; // Paul
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
         printf("Our alpha is %8.4f\n", alpha);
         printf ("         n     C**bar        +/-        t       t*\n");
      }

   while (alphafound != 1){
      A2bar = XAbar = X2bar = 0;
      for (int n = 1; n <= 100000; ++n) {
         // Initialize the stock price.
         S = S0;

         // Initialize the Brownian path.
         B = 0;

         // Initialize time.
         t = 0;

         for (int i = 1; i <= N; ++i) {

               // Advance the path.
               U = MTUniform(0);
               Z = PsiInv(U); // Standard normal via inverst transform.
               B += sqrtdt * Z; // standard brownian path 
               Btilde = -B; // symmetric distribution

               // Advance time by one period.
               t += dt;

               // Compute the next stock price.
               S = S0 * exp(mu * t + sigma * B);            // standard stock motion
               Stilde = S0 * exp(mu * t + sigma * Btilde);  // antithetic stock process

               //// Compute the value of A (control variable) 
               A = B * B - t;

               // Determine call payoff
               C = max(S - K);
               Ctilde = max(Stilde - K);

               // Determine C*
               X = exp(-r * t) * ((C + Ctilde) / 2);
               A2bar = ((n - 1) * A2bar + A * A) / n;
               Xbar = ((n - 1) * Xbar + X) / n;
               X2bar = ((n - 1) * X2bar + X * X) / n;
               XAbar = ((n - 1) * XAbar + X * A) / n;
         }
      }
      alpha = -XAbar / A2bar;
      CorA_Cstar = XAbar / ( sqrt(X2bar - Xbar*Xbar) * sqrt(A2bar) );

      alphafound = 1;
   }

      // Start the clock to time the computations.
      timestart = clock ();

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


            A = B*B - T;

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
               //printf("   This alpha calculated is %8.3f \n", alpha);
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
               part6alpha = alpha;
            }
            else if (part == 7) {
               printf ("%10.0f   %8.4f   %8.6f %8.3f %8.3f\n", n, Cstarstarbar, error, elapsed_time, t_star);
               part7 = Cstarstarbar;
               part7alpha = alpha;
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

   printf("\nPart 6:\n=======\n   The price a call option is %3.4f \n", part6);
   printf("   This was significantly faster than Part 3, running in %8.3f seconds\n", part6time);
   printf("   This alpha calculated is %8.3f \n", part6alpha);
   printf("   This sample Cov(A,C*) was calculated to be %8.3f \n", XAbar);
   printf("   This sample Cor(A,C*) was calculated to be %8.3f \n", CorA_Cstar);

   printf("\nPart 7:\n=======\n   The value of the look-back option is %2.5f\n", part7);
   printf("   This alpha calculated is %8.3f \n", part7alpha);



   Exit ();
   
}



