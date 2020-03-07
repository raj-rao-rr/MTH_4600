#include <iostream>
#include "Functions.h"

double max(double val) {
    if (val > 0.0) {
        return val;
    }
    return 0.0;
}

int main() {

   double r, t, T, mu, sigma, alpha, dt, sqrtdt, S, Stilde, S0, V, Vbar, V2bar,
          elapsed_time, t_star, stdhatV, error, epsilon, n, Discount_factor,
          U, Z, A, B, Btilde, C, Ctilde, BSprice, K=110.0;
   int i, N, done, test;

   // appropriate value for alpha 
   alpha = 0.1;

   // Time to expiration.
   T = 0.5;

   // Number of stock price periods.
   N = 50;

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

   // Specify the 95% error tolerance.
   epsilon = 0.005;

   // Start the clock to time the computations.
   Time ();

   // Seed the RNG.
   MTUniform (1);

   // Print column headings for output to execution window.
   printf ("         n       Vbar        +/-        t       t*\n");

   // Initialize certain values.
   V2bar = Vbar = done = n = test = 0;

   // Begin the simulation loop.
   while (!done) {

      // Initialize the stock price.
      S = S0;

      // Initialize the Brownian path.
      B = 0;

      // Initialize time.
      t = 0;

      // Simulate a stock price path.  Go forward period-by-period computing
      //   the next stock price.
      for (i = 1; i <= N; i++) {

         // Advance the path.
         U = MTUniform (0);
         Z = PsiInv (U); // Standard normal via inverst transform.
         B += sqrtdt * Z; // standard brownian path 
         Btilde = -B; // symmetric distribution

         // Advance time by one period.
         t += dt;

         // Compute the next stock price.
         S = S0 * exp (mu*t + sigma*B);  // standard stock motion
         Stilde = S0 * exp(mu * t + sigma*Btilde); // antithetic stock process


      }
      // S now is the value of the stock at time T for the simulated price path.

      // Actual Price of the Option (problem 5)
      BSprice = BlackScholes(T, S0, K, sigma, r);

      // Discount back to time 0.
      C = max(S-K);
      Ctilde = max(Stilde - K);
      A = B*B - T;
      V = (Discount_factor * (C + Ctilde)/2) + (alpha * A);

      // Update the simulation counter and the sample moments of V at time T.
      n ++;
      Vbar  = ((n-1)*Vbar + V) / n;
      V2bar = ((n-1)*V2bar + V*V) / n;

      // Update the error condition test counter.
      test ++;

      // Test the error condition periodically.
      if (test == 100000) {

         // Estimate the standard deviation of S.
         stdhatV  = sqrt (V2bar - Vbar*Vbar);

          // Estimate the error of the Vbar estimator.
         error = 1.96 * stdhatV / sqrt(n);

         // Compute the elapsed time and estimate the time of completion.
         elapsed_time = Time ();
         t_star = elapsed_time * pow (error / epsilon, 2);

         // Report.
         printf ("%10.0f   %8.4f   %8.6f %8.3f %8.3f\n", n, Vbar, error, elapsed_time, t_star);
        
         // Reset the "test" counter and see if error tolerance is met.
         test = 0;
         if (error < epsilon) {
            done = 1;
         }

      }
   }
   printf("The exact value of the option contract is %8.4f\n", BSprice);
   // printf("Correlation between C* and A is %8.4f\n", Correlation());
   Exit ();

}



