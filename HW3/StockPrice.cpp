#include <iostream>
#include "Functions.h"

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

   double r, t, T, mu, sigma, alpha, dt, sqrtdt, S, Stilde, S0, V, Vbar, V2bar,
          elapsed_time, t_star, stdhatV, error, epsilon, n, Discount_factor,
          U, Z, A, B, Btilde, C, Ctilde, BSprice, A2bar, XAbar, K, Smax, SmaxT;
   int i, N, done, test;

   // Strike Price
   K = 110.0;

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

      // Maximum stock price 
      Smax = SmaxT = 0;

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
         U = MTUniform (0);
         Z = PsiInv (U); // Standard normal via inverst transform.
         B += sqrtdt * Z; // standard brownian path 
         Btilde = -B; // symmetric distribution

         // Advance time by one period.
         t += dt;

         // Compute the next stock price.
         S = S0 * exp (mu*t + sigma*B);            // standard stock motion
         Stilde = S0 * exp(mu * t + sigma*Btilde); // antithetic stock process

         //// Compute the value of A
         A = B * B - T;

         // Update the required sample moments.
         A2bar += A * A;
         XAbar += B * A;

         // Keep track of the maximum stock price 
         if (S > Smax) {
             Smax = S;
         }
         if (Stilde > SmaxT) {
             SmaxT = Stilde;
         }

      }
      // S now is the value of the stock at time T for the simulated price path.

      // Actual Price of the Option (problem 5)
      BSprice = BlackScholes(T, S0, K, sigma, r);

      // Determine call payoff
      C = max(Smax - K);
      Ctilde = max(SmaxT - K);

      // Appropriate value for alpha 
      alpha = -XAbar/ A2bar;

      // Sample mean statistics being averaged 
      V = Discount_factor * ((C + Ctilde)/2) + (alpha * A);

      // Update the simulation counter and the sample moments of V at time T.
      ++ n;
      Vbar  = ((n-1)*Vbar + V) / n;
      V2bar = ((n-1)*V2bar + V*V) / n;

      // Update the error condition test counter.
      ++ test;

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
   printf("Our alpha term is %8.4f\n", alpha);
   Exit ();

}



