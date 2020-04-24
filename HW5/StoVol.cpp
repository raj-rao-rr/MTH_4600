////////////////////////////////////////////////////////////////////////////////
// This program generates 10 stock price paths using a GARCH(1,1)
//    stochastic volatility model. Results are printed to the screen. The
//    current (stochastic) and long-run volatility are 30% annually.
////////////////////////////////////////////////////////////////////////////////

#include "Functions.h"

int main () {

   int i, n, N;
   double s_start, s, T, r, sigma, sigma2,  sigma2_start, mu, alpha, beta,
          gamma, dt, N01, R, U, V2;

   // alpha, beta, and gamma parameters as computed in the XOM example.
   alpha = 0.045; beta = 0.083; gamma = 0.872;

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

   // Generate 10 illustrative paths.
   for (n = 1; n <= 10; n++) {

      // Initialize the stock price and the volatility to their starting (time 0)
      //    values.
      s = s_start;
      sigma2 = sigma2_start;

      // Now march forward in time day by day.
      for (i = 1; i <= N; i++) {

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
         printf ("%7.3f  %7.3f  %7.3f\n", i*dt, s, 100.0 * sqrt(sigma2 * 252.0));

      }

      Pause ();

   } // Next path.

}


