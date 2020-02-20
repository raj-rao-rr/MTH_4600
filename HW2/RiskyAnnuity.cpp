/*******************************************************************************
This code values an annuity that makes "promised" payments of $100 per month for
30 years.  The annuity is valued for the following 6 types of credits:
    - 0 is risk free, never defaults;
    - 1 is AAA, time to default (survival time) is Exponential (200 years);
    - 2 is AA,  survival time is Exponential (100);
    - 3 is A,   survival time is Exponential (66 2/3);
    - 4 is BBB, survival time is Exponential (50);
    - 5 is "Speculative", survival time is Exponential (33 1/3);
The program calculates the present value (at 3% compounded continuously) of
the stream of expected cash flow, as discussed in class.
*******************************************************************************/

double YTM (double, double, int);

#include "Functions.h"



int main() {

   int credit, m, maturity;
   double lambda[] = {1000000.0, 200.0, 100.0, 66.666667, 50.0, 33.333333};
   double rate = 0.03, s, pv, pv0, cashflow;

   // Thirty year maturity in months.
   maturity = 360;

   // Promised monthly cash flow is $100.
   cashflow = 100.0;

   // Do the present value calculation for each credit.
   for (credit = 0; credit <= 5; credit ++) {

      // Initialize the present value to 0.
      pv = 0.0;

      // Loop through the monthly payments to maturity.
      for (m = 1; m <= maturity; m ++) {

         // Calculate the probability of survival to month number "month".
         // Lambda is in years.
         s = exp (-m / (12.0*lambda[credit]));

         // Augment the pv by this month's discounted expected cash flow.
         pv += exp (-rate * m / 12.0) * s * cashflow;

      }

      // Record the pv for riskless cash flow.
      if (credit == 0) {
         pv0 = pv;
      }

      // Output the results to the screen.
      printf ("For credit number %d, the present value of the expected cash flow is %8.2f.\n", credit, pv);
      printf ("As a percent of the riskless PV, this is %8.2f.\n", 100.0 * pv/pv0);
      printf ("The yield-to-maturity is %8.2f percent.\n\n\n", YTM (pv, cashflow, maturity));


   }

   Pause ();

}

////////////////////////////////////////////////////////////////////////////////
// This function calculates the yield-to-maturity (ytm) of the promised cash flow
//   of the annuity where the calculated value of this promised cash flow
//   is "pv". The ytm is calculated as a continuously compounding annual rate.
////////////////////////////////////////////////////////////////////////////////
double YTM (double pv, double cashflow, int maturity) {

   int m;
   double pv0, ytm = 0.0, step = 0.01;

   // Calculate the ytm to within 0.000001, i.e., 0.0001% or 0.01 bps
   while (step > 0.000001) {

      // Keep increasing the trial ytm until the resulting present value (pv0)
      //  is too small (less than pv).
      while (1) {

         // Calculate the present value with a discount rate of "ytm".
         pv0 = 0;
         for (m = 1; m <= maturity; m++) {
            pv0 += cashflow * exp (-ytm * m / 12.0);
         }

         // Is ytm now too big? If not, increase it further by the amount "step".
         if (pv0 > pv) {
            ytm += step;
         }

         // If so, reduce it by "step" and break out of the "while (1)" loop.
         // Reduce the step size for the next iteration.
         else {
            ytm -= step;
            step /= 10.0;
            break;
         }

      }

   }

   // Return the yield-to-maturity as a percent.
   return 100.0 * ytm;

}               


