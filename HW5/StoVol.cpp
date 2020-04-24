////////////////////////////////////////////////////////////////////////////////
// This program generates 10 stock price paths using a GARCH(1,1)
//    stochastic volatility model. Results are printed to the screen. The
//    current (stochastic) and long-run volatility are 30% annually.
////////////////////////////////////////////////////////////////////////////////

#include "Functions.h"

// Used for calculating the max payoff of a simple vanilla call
double Max(double &val) {
    if (val > 0.0) { return val; }
    return 0.0;
}

int main () {

   int i, n, N, sims;
   double s_start, s, T, r, sigma, sigma2,  sigma2_start, mu, alpha, beta,
          gamma, dt, N01, R, U, V2;
   double val1, val2, avg1 = 0;

   double strikes[11] = { 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0 };
   double averages[11] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

   // alpha, beta, and gamma parameters as computed in the XOM example.
   alpha = 0.045; beta = 0.083; gamma = 0.872;

   // number of simulations to run
   sims = 100000;

   // Time to expiration.
   T = 0.5;

   // Risk-free interest rate.
   r = 0.05;

   // Number of days to expiration (252 days per year), and one day in years.
   N = 252 * T;
   dt = 1.0/252.0;

   // Convert "r" to a daily risk-free rate.
   r *= dt;

   // Current stock price.
   s_start = 100;

   // Current stock price variance.
   sigma2_start = 0.35 * 0.35;

   // Current daily stock price variance.
   sigma2_start *= dt;

   // Annual long-term variance.
   V2 = 0.30 * 0.30;

   // Daily long-term variance.
   V2 *= dt;

   // Seed the RNG.
   MTUniform (1);

   // Generate Monte Carlo paths for Stochastic volatility
   for (n = 1; n <= sims; ++n) {

      // Initialize the stock price and the volatility to their starting (time 0)
      //    values.
      s = s_start;
      sigma2 = sigma2_start;

      // Now march forward in time day by day.
      for (i = 1; i <= N; ++i) {

         // Compute the drift parameter that corresponds to the volatility for
         //    this period.
         mu = r - 0.5*sigma2;

         // Compute the stock price at day "i"...

         // First get a standard normal N01.
         U = MTUniform (0);
         N01 = PsiInv (U);

         // Apply current volatility.
         sigma = sqrt(sigma2);
         R = sigma * N01;

         // Update the stock price.
         s *= exp (mu + R);

         // Update the stochastic volatility according to the GARCH(1,1) model.
         sigma2 = alpha * V2  +  beta * R*R +  gamma * sigma2;

         // Print the current stock price and annualized volatility to
         //    the screen.
         // printf ("%7.3f  %7.3f  %7.3f\n", i*dt, s, 100.0 * sqrt(sigma2 * 252.0));

      }

      // Problem 1 //
      
      // check the payoff with stock at terminal time t
      val1 = s - 0;
      avg1 += Max(val1);

      // Problem 2 //
      for (int j = 0; j < 11; ++j) {
          val2 = s - strikes[j];
          averages[j] += Max(val2);
      }

   } 

   // Outputs value for problem 1
   printf("Problem 1\n");
   printf("Our call option with strike K = 0 is valued at %8.4f\n", avg1 / sims);
   printf("--------------------------------------------------------------\n");

   // Outputs value for problem 2
   printf("Problem 2\n");
   for (int t = 0; t < 11; ++t) {
       printf("Our call option with strike K = %8.2f is valued at %8.4f\n", strikes[t], averages[t] / sims);
   }
   printf("--------------------------------------------------------------\n");

   Exit();
}


