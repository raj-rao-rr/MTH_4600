////////////////////////////////////////////////////////////////////////////////
// Starter code for Homework 6.
////////////////////////////////////////////////////////////////////////////////


// These functions are found below.
void GetData ();
void Report ();

// These global variables are used to reduce function arguments.
double *R, *R2, *sigma;

#include "Functions.h"

int main () {

   double *s, *p, *ptilde, *g;
   double sigma0, alpha;
   double *AllocateMemory();

   // Allocate necessary memory.
   s = AllocateMemory ();
   p = AllocateMemory ();
   ptilde = AllocateMemory ();
   g = AllocateMemory ();


   // Read in the time series data.
   GetData ();

   // Implement the hidden Markov chain model between here...

   // assigning the free parameters 
   sigma0 = 0.25;
   alpha = 0.03;

   for (int i = -1258; i <= 1258; i += 2) {
       // the daily volatility when the walk is in state k
       s[1259 + i] = sigma0 * exp(alpha * i);

       // 
   }

   for (int i = 1; i <= 1258; ++i) {
       printf("%8.4f\n", R[i]);
   }



   // ... and here.
   // !!! When you have done this, get read "****" below in the Report() function.

   // Create TeX files for viewing results.
   Report ();
   Exit();
}

////////////////////////////////////////////////////////////////////////////////
// Allocate memory for an array with indices from -1260 to +1260.
// This is a little more than needed.
////////////////////////////////////////////////////////////////////////////////

double *AllocateMemory () {

   double *x;

   x = (double *) calloc (2*1260 + 1, sizeof (double));

   x += 1260;

   return x;

}

////////////////////////////////////////////////////////////////////////////////
// Read in a daily time series of stock price returns R[t] 1 <= t <= 1258.
////////////////////////////////////////////////////////////////////////////////
void GetData () {

   int t;
   char input[100];
   FILE *fp;

   fp = fopen ("XOM5YrsDaily.txt", "r");

   // Read in the file description.
   fgets (input, 99, fp);

   // Allocate memory for the data, initialized to 0.
   R     = List (1258);
   R2    = List (1258);
   sigma = List (1258);

   // Now read in the data.
   for (t = 1; t <= 1258; t++) {
      fgets (input, 99, fp);
      sscanf (input, "%lf", R+t);
      R2[t] = R[t] * R[t];
   }

   fclose (fp);

   return;

}

////////////////////////////////////////////////////////////////////////////////
// Generate some output files.
////////////////////////////////////////////////////////////////////////////////
void Report () {

   int t;
   double annualizedVol, Z;
   FILE *fp1, *fp2;


   fp1 = fopen ("HistEmpVols.txt", "w");
   fp2 = fopen ("StandardizedXOM.txt", "w");

   // Start at day 51.
   for (t = 51; t <= 1258; t++) {

      annualizedVol = sqrt (252.0) * sigma[t];
      fprintf (fp1, "%4d %10.2f\n", t, annualizedVol);

      // Standardize the data and add to a histogram.
      Z = PsiInv(MTUniform(0));               // **** Get rid of this line,
      //Z = R[t] / sigma[t];                  // and un-comment this line.
      NormalHistogram (Z, 20, 0);
      fprintf (fp2, "%4d %10.2f\n", t, Z);

   }

   fclose (fp1);
   fclose (fp2);

   // Create histogram TeX file.
   NormalHistogram (0, 20, 1);

   return;

}






