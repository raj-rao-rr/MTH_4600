#include <iostream>
#include "Functions.h"
#include <string>
using namespace std; 

// Used for calculating the max payoff of a simple vanilla call.
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
   // Declare all variables.
   double r, t, T, mu, sigma, alpha, alphafound, dt, sqrtdt, S, Stilde, S0, V, Vbar, V2bar,
          elapsed_time, timestart, t_star, stdhat, error, epsilon, n, Discount_factor,
          U, Z, A, B, Btilde, C, Ctilde, Cstar, Cstarbar, Cstar2bar,Cstarstar, Cstarstarbar, Cstarstar2bar,
          BSprice, A2, XA, K, Smax, SmaxT, part1, part1time, part2, part2time, part3, part3time, part4, part5,
          part6, part6time, part6alpha, part7, part7alpha, part8, X, A2bar, Xbar, XAbar, X2bar, CorA_Cstar,
          part6Cov, Stau, shout_avg, part8Vbar;
   int i, N, done, test, part;

   // Similar to the variable "done," alphafound indicates when alpha has been calculated
   // so that it is only calculated once.
   alphafound = 0;

   // This c++ file goes through every part of homework 3.
   for(part = 1; part <= 8; ++part){
      printf("\n\nPart %d:\n", part);
      
      // For testing purposes, some of the more computationally and time heavy parts can be skipped over.
      if (part == 0 || part == 0 || part == 0 || part == 0 || part == 0) { 
         printf("Part %d omitted for testing purposes\n", part);
         continue;
      }

      // Strike Price of the option.
      K = 110.0;

      // Time to expiration.
      T = 0.5;

      // Number of stock price periods. The standard for this homework is 50 periods.
      // If Part 4 is being computed, the number of periods is switched to 1.
      N = 50;
      if (part == 4) {
         N = 1;
      }

      // Time increment per period.
      dt = T / N;

      // Risk-free interest rate.
      r = 0.05;

      // Compute these oft-used values once and for all.
      Discount_factor = exp(-r*T);
      sqrtdt = sqrt(dt);

      // Stock price volatility.
      sigma = .30;

      // Drift term.
      mu = r - sigma*sigma/2.0;

      // Initial stock price.
      S0 = 100.0;

      // If Part 5 is being computed, only the Black Scholes formula needs to be used
      // so the Black Scholes function from Functions.h is used and this iteration of
      // the "part loop" is ended with "continue."
      if (part == 5) { 
         // Actual Price of the Option calculated in part 5.
         part5 = BlackScholes(T, S0, K, sigma, r);
         printf("   This uses the Black-Scholes call option pricing formula.");
         continue;
      }

      // Specify the 95% error tolerance.
      epsilon = 0.005;

      // Seed the RNG.
      MTUniform (1);

      // Print column headings for output to execution window. These column
      // headings change depending on the part of the homework being computed.
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

   // This while loop is used to determine if the value alpha has been caluclated or not.
   // If alphafound is 0, then we enter the loop, otherwise, it is skipped.
   while (alphafound != 1){

      // Initialized certain variables.
      A2bar = XAbar = X2bar = Xbar = 0;

      for (int n = 1; n <= 100000; ++n) {
         // Initialize the stock price.
         S = S0;

         // Initialize the Brownian path, and the time.
         B = t = 0;

         for (int i = 1; i <= N; ++i) {
               // Construct the Brownian path.
               U = MTUniform(0);
               Z = PsiInv(U);    // Standard normal via inverst transform.
               B += sqrtdt * Z;  // standard brownian path.
               Btilde = -B;      // symmetric distribution.

               // Advance time by one period.
               t += dt;

               // Compute the next stock price.
               S = S0 * exp(mu * t + sigma * B);            // Standard stock motion.
               Stilde = S0 * exp(mu * t + sigma * Btilde);  // Antithetic stock process.

               // Compute the value of A which is the control variable.
               A = B * B - t;

               // Determine call payoff.
               C = max(S - K);
               Ctilde = max(Stilde - K);
               
               // Calculate C* (here called X).
               X = exp(-r * t) * ((C + Ctilde) / 2);

               // Calculate the Cov(A,C*) and the var(A) needed to find alpha.
               A2bar = ((n - 1) * A2bar + A * A) / n;
               XAbar = ((n - 1) * XAbar + X * A) / n;

               // Calculate Xbar and X2bar to be able to calulcate the stdev of X
               // for the Cor(A,C*) calculation.
               Xbar = ((n - 1) * Xbar + X) / n;
               X2bar = ((n - 1) * X2bar + X * X) / n;
         }
      }
      // alpha and Cor(A,C*) are calculated once and for all.
      // There is no need to calculate Abar, because it is known to be 0.
      alpha = - XAbar / A2bar;
      part6Cov = XAbar;
      CorA_Cstar = XAbar / ( sqrt(X2bar - Xbar*Xbar) * sqrt(A2bar) );

      // alphafound is assigned 1, so this loop will no longer be repeated.
      alphafound = 1;
   }

      // Start the clock to time the computations.
      timestart = clock ();

      // Initialize certain values.
      V2bar = Vbar = A2 = XA = done = n = t = test = 0;
      Cstarstar2bar = Cstarstar2bar = Cstar2bar = Cstarbar = 0; // Paul

      // Begin the simulation loop.
      while (!done) {

         // Initialize the stock price with the intial stock price.
         S = S0;

         // Initialize max stock prices, Brownian path, and time.
         Smax = SmaxT = B = t = 0;

         // Simulate a stock price path.  Go forward period-by-period computing
         // the next stock price.
         for (i = 1; i <= N; ++i) {
            
            // Advance time by one period.
            t += dt;

            // Advance the path.
            U = MTUniform (0);
            Z = PsiInv (U); // Standard normal via inverse transform.
            B += sqrtdt * Z; // Standard Brownian path 
            S = S0 * exp (mu*t + sigma*B); // Compute the next stock price.
            
            // Only Part 1 does not use the antithetic variance reduction strategy,
            // so if the Part is not 1, then calculate Btilde (which has the same 
            // distribution as B), and find the antithetic stock price path Stilde.
            if (part != 1){
               Btilde = -B; // symmetric distribution
               Stilde = S0 * exp(mu * t + sigma*Btilde); // antithetic stock process
            }

            // Part 7 deals with the "look-back” option, therefore the maximum value
            // of S along its path must be recorded. This section does that as follows:
            // if a new maximum is found, replace the old maximum, otherwise do nothing.
            if (part == 7) {
               if (S > Smax) {
                  Smax = S;
               }
               if (Stilde > SmaxT) {
                  SmaxT = Stilde;
               }
            }
         }
         
         // S now is the value of the stock at time T for the simulated price path.

         // Discount the simulated price of S at time T back to the present. 
         if (part == 1) {
            V = Discount_factor * S;
         }

         // Discount the simulated price of S at time T (this time obtained using the 
         // antithetic variance reduction technique) back to the present.
         else if (part == 2) {
            V = Discount_factor * (S + Stilde) / 2.0 ;
         }

         // Parts 3 and 4 require the price of a call option to be calculated.
         // Recall that for Part 4, only on period is used, rather than 50.
         else if (part == 3 || part == 4) { 

            // Determine call payoff at time T (antithetic method is used here as well).
            C = max(S - K);
            Ctilde = max(Stilde - K);

            // Calculate C*, which is the payoff discounted to the present.
            Cstar = Discount_factor * (C + Ctilde)/2;
         }

         // Parts 6 and 7 also require the price of a call option to be calculated.
         // Diffent from Parts 3 and 4, a control varible A is used to reduce the variance
         // and thus shorten the run-time.
         else if (part == 6 || part == 7) {

            if (part == 6){
               // Part 6 uses the last stock price observed, just like in Parts 3 and 4.
               C = max(S - K);
               Ctilde = max(Stilde - K);
            } 

            else if (part == 7){
               // Part 7 uses the maximum stock price observed along the stock price path,
               // here represented by Smax, and SmaxT.
               C = max(Smax - K);
               Ctilde = max(SmaxT - K);
            }

            // The control variable A is computed. It has an expected value of 0.
            A = B*B - T;

            // Here, C** is calculated. It is the payoff of the option at time T 
            // (calculated with the same antithetic method as before), with an additional
            // control varibal added to further reduce the variance and run-time.
            // Alpha was previously  calculated at the onset of the code.
            Cstarstar = Discount_factor * ((C + Ctilde)/2) + (alpha * A);
         }

         // Update the simulation counter and the sample moments of V at time T.
         ++ n;

         // Different parts ask to calculate/estimate prices for different instruments.
         // We must later use these prices and their variance to determine the error condition.
         // Parts 1 and 2 ask for Vbar, the average discounted terminal stock price.
         if (part == 1 || part == 2) {
            Vbar  = ((n-1)*Vbar + V) / n;
            V2bar = ((n-1)*V2bar + V*V) / n;
         }
         // Parts 3 and 4 ask for C*, the average discounted call option payoff
         // with antithetic reduction applied.
         else if (part == 3 || part == 4) {
            Cstarbar  = ((n-1)*Cstarbar + Cstar) / n;
            Cstar2bar = ((n-1)*Cstar2bar + Cstar*Cstar) / n;
         } 
         // Part 6 asks for C**, the average discounted call option payoff.
         // Part 7 asks for C**, the average discounted look-back option payoff.
         // Both use antithetic reduction and control variable A.
         else if (part == 6 || part == 7) {
            Cstarstarbar  = ((n-1)*Cstarstarbar + Cstarstar) / n;
            Cstarstar2bar = ((n-1)*Cstarstar2bar + Cstarstar*Cstarstar) / n;
         }

         // Update the error condition test counter.
         ++ test;

         // Test the error condition periodically.
         if (test == 100000) {

            // In Parts 1 and 2, the error condition is based on Vbar.
            if (part == 1 || part == 2){
               // Estimate the standard deviation of Vbar.
               stdhat  = sqrt (V2bar - Vbar*Vbar);
            }

            // In Parts 3 and 4, the error condition is based on C*bar.
            else if (part == 3 || part == 4){
               // Estimate the standard deviation of C*bar.
               stdhat  = sqrt (Cstar2bar - Cstarbar*Cstarbar);
            }

            // In Parts 6 and 7, the error condition is based on C**bar.
            else if (part == 6 || part == 7){
               // Estimate the standard deviation of C**bar.
               stdhat  = sqrt (Cstarstar2bar - Cstarstarbar*Cstarstarbar);
            }

            // Estimate the error of the Vbar/C*bar/C**bar estimator.
            error = 1.96 * stdhat / sqrt(n);

            // Compute the elapsed time and estimate the time of completion.
            elapsed_time = (double) (clock() - timestart) / CLOCKS_PER_SEC;
            t_star = elapsed_time * pow (error / epsilon, 2);

            // Reset the "test" counter and see if error tolerance is met.
            test = 0;
            if (error < epsilon) {
               done = 1;
            }

            // This section 1) Reports and 2) Saves the important outputs to then print them all at the end.
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

      if (part == 8){

         // Average of Shout Option Pricing (varied tau's)
         shout_avg = 0;

         // Start the clock to time the computations.
         Time ();

         // Print column headings for output to execution window.
         printf ("         n       Vbar        +/-        t       t*\n");

         // Initialize certain values.
         V2bar = Vbar = done = n = test = 0;

         // succesively itterate for values of tau
         for (int tau = 0; tau <= 0; ++tau) {
            // printf("At tau = %8.4d \n", tau);

            // Begin the simulation loop.
            while (!done) {

               // Initialize the stock price.
               S = S0;

               // Initialize the Brownian path.
               B = 0;

               // Initialize time.
               t = 0;

               // A: Initial phase: estimate Var A, Cov (X,A).
               A2bar = XAbar = 0;

               // Simulate a stock price path.  Go forward period-by-period computing
               //   the next stock price.
               for (i = 1; i <= N; ++i) {

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

                     // Keep track of S-tau
                     if (tau == 0) {
                        Stau = S0;
                     }
                     else {
                        if (tau == i) {
                           Stau = S;
                        }
                     }
               }
               // S now is the value of the stock at time T for the simulated price path.

               // Determine call payoff
               C = max(S - K);
               Ctilde = max(Stilde - K);

               // Determine C* (anthetic reduction call price) 
               X = exp(-r * t) * ((C + Ctilde) / 2);

               // Compute the value of A (control variable) 
               A = B * B - T;

               // Sample mean statistics being averaged 
               V = X + (alpha * A);

               // Update the simulation counter and the sample moments of V at time T.
               ++n;
               Vbar = ((n - 1) * Vbar + V) / n;
               V2bar = ((n - 1) * V2bar + V * V) / n;

               // Update the error condition test counter.
               ++test;

               // Test the error condition periodically.
               if (test == 100000) {

                     // Estimate the standard deviation of S.
                     stdhat = sqrt(V2bar - Vbar * Vbar);

                     // Estimate the error of the Vbar estimator.
                     error = 1.96 * stdhat / sqrt(n);

                     // Compute the elapsed time and estimate the time of completion.
                     elapsed_time = Time();
                     t_star = elapsed_time * pow(error / epsilon, 2);

                     // Report.
                     printf("%10.0f   %8.4f   %8.6f %8.3f %8.3f\n", n, Vbar, error, elapsed_time, t_star);
                     // Reset the "test" counter and see if error tolerance is met.
                     test = 0;
                     if (error < epsilon) {
                        done = 1;
                     }

               }
            }
            part8Vbar = Vbar;
            shout_avg += Vbar;
         
         }

         shout_avg /= N;
      }
   }

   // All the outputs saved from the iterations of the "part loop" that have been saved
   // are now pulled out to print the answers to every part.
   printf("\nPart 1:\n=======\n   Vbar is %3.4f and therefore agrees with S0\n", part1);
   printf("   This part ran in %8.3f seconds\n", part1time);

   printf("\nPart 2:\n=======\n   Similar to Part 1, Vbar is %3.4f and therefore also agrees with S0\n", part2);
   printf("   This was significantly faster than Part 1, running in %8.3f seconds\n", part2time);

   printf("\nPart 3:\n=======\n   Modifying the code of Part 2, the price a call option is %2.5f\n", part3);
   printf("   This part ran in %8.3f seconds\n", part3time);

   printf("\nPart 4:\n=======\n   Redoing Part 3 with N=1 instead of 50, the price a call option is %2.5f\n", part4);
   printf("   This is the same answer as in Part 3\n");

   printf("\nPart 5:\n=======\n   Using the Black-Scholes formula, the exact value of the option is %2.5f\n", part5);

   printf("\nPart 6:\n=======\n   The price a call option is %3.4f \n", part6);
   printf("   This was significantly faster than Part 3, running in %8.3f seconds\n", part6time);
   printf("   The alpha calculated was calculated to be %8.3f \n", part6alpha);
   printf("   The sample Cov(A,C*) was calculated to be %8.3f \n", part6Cov);
   printf("   The sample Cor(A,C*) was calculated to be %8.3f \n", CorA_Cstar);

   printf("\nPart 7:\n=======\n   The value of the look-back option is %2.5f\n", part7);
   printf("   The alpha was calculated to be %8.3f \n", part7alpha);
   
   printf("\nPart 8:\n=======\n   The shout option is priced at %8.4f \n", shout_avg);
   printf("   Our Option is priced at %8.4f \n", part8Vbar);
   printf("---------------------------------------------------------------\n");


   // We're done
   Exit ();
   
}



