////////////////////////////////////////////////////////////////////////////////
// This program calculates the value of call options via Black-Scholes.  The
//    time to expiration is fixed at 1/2 year; the risk-free rate at 5%; and
//    the current value of the underlying at $100.
////////////////////////////////////////////////////////////////////////////////

#include "Functions.h"

int main () {

   double s0, T, r, sigma, c, k;

   printf ("This program calculates via Black-Scholes the value of\n");
   printf ("a call option expiring in 0.5 years when the underlying\n");
   printf ("stock price is 100 and the risk-free rate is 5%%.\n\n");


   // Time to expiration.
   T = 0.5;

   // Current stock price.
   s0 = 100;

   // Risk-free interest rate.
   r = 0.05;

again:

   // Get the stock price volatility.
   sigma = GetDouble ("What is the volatility in percent?... ");
   sigma /= 100.0;

   // Get the option's strike price.
   k = GetDouble ("What is the strike?...                ");

   c = BlackScholes (T, s0, k, sigma, r);
   printf ("The option's value is:               %8.4f\n\n", c);


   goto again;

}




