////////////////////////////////////////////////////////////////////////////////
// This program calculates the implied Black-Scholes volatility of call options.
// The time to expiration is fixed at 1/2 year; the risk-free rate at 5%; and
//    the current value of the underlying at $100.
////////////////////////////////////////////////////////////////////////////////

#include "Functions.h"

int main () {

   double s0, T, r, sigma, c, k;

   printf ("This program calculates the Black-Scholes implied volatility\n");
   printf ("for a call option expiring in 0.5 years when the underlying\n");
   printf ("stock price is 100 and the risk-free rate is 5%%.\n\n");

   // Time to expiration.
   T = 0.5;

   // Current stock price.
   s0 = 100;

   // Risk-free interest rate.
   r = 0.05;

again:

   // Get the option value.
   c = GetDouble ("What is the option price?... ");

   // Get the option's strike price.
   k = GetDouble ("What is the strike?...       ");

   sigma = ImpliedVol (T, s0, k, r, c);
   printf ("The implied annual volatility in percent is:  %8.4f\n\n", 100.0*sigma);


   goto again;

}



