#include <iostream>
#include <list>
#include "Functions.h"

// MTUniform is the Mersenne twister, MWCUniform is the multiply with carry
//    generator that we discussed in class.


double MWCUniform (unsigned int);
unsigned int Temper (unsigned int);

int main() {

   std::cout << "\n\n\n\n\n\n";
   // Insert your code between here...
   int M1, M2, M3, Triplet, Trials;
   double p, q, mu, sigma, Z;
   std::list<double> listOfFloats;
   int X[64000]= { };

   Trials = 100000000;

	// Create Triplets with the bin that they're in, then update the counter
	for (int i = 0; i < Trials; i++){
      M1 = int(MWCUniform(0) / 0.025) * 40 * 40;
      M2 = int(MWCUniform(0) / 0.025) * 40;
      M3 = int(MWCUniform(0) / 0.025);
      Triplet = M1 + M2 + M3;
      ++X[Triplet];
		listOfFloats.push_back(Triplet);
   }

   // Compute the mean and standard deviation of the X[]'s.
   p = 1.0/64000.0;
   q = 1.0 - p;
   mu = Trials * p;
   sigma = sqrt(Trials * p * q);

   for (int i = 0; i < 64000; i++) {
      Z = (X[i] - mu) / sigma;
      NormalHistogram (Z, 40, 0);
   }

   // Create output files.
   NormalHistogram (0, 40, 1);


   // ... and here.
   std::cout << "\n\n\n\n\n\n";

   Exit ();

}



////////////////////////////////////////////////////////////////////////////////
// The following code implements the Multiply-With-Carry RNG.
////////////////////////////////////////////////////////////////////////////////
double MWCUniform (unsigned int seed) {

   // Static variables retain their values between function calls.
   static unsigned int a1, a0, n1, n0, c1, c0, initialized=0;

   unsigned int w1, w0, x1, x0, y1, y0, z1, z0, carry, N;

   // This function is called only be MWCUniform and follows below.
   void Split (unsigned int, unsigned int *, unsigned int *);

   // Seed the MWC function when it is passed a non-zero seed or it is not yet
   //    initialized.
   // Represent a, n, and c as with pairs of 16-bit numbers:
   //  a = a1 * 2^16 + a0 (the multiplier)
   //  n = n1 * 2^16 + n0 (initially the seed)
   //  c = c1 * 2^16 + c0 (the increment, initially 1)
   if (seed || !initialized) {
      seed = (seed ? seed : 1);
      Split (4294967118, &a1, &a0);
      Split (seed, &n1, &n0);
      c0 = 1;                                 // first 16 bits of "c"
      c1 = 0;                                 // second 16 bits
      initialized = 1;
   }

   // Generate the next number in the sequence.

   // Compute w1, w0, where a0 * n0 = w1 * 2^16 + w0.
   Split (a0 * n0, &w1, &w0);

   // Compute x1, x0, where a1 * n0 = x1 * 2^16 + x0.
   Split (a1 * n0, &x1, &x0);

   // Compute y1, y0, where a0 * n1 = y1 * 2^16 + y0.
   Split (a0 * n1, &y1, &y0);

   // Compute z1, z0, where a1 * n1 = z1 * 2^16 + z0.
   Split (a1 * n1, &z1, &z0);

   // Compute next n0, n1, c0, c1, where an + c = c1 * 2^48 + c0 * 2^32 + n1 * 2^16 + n0.

   // First compute n0 and whatever is carried.
   Split (w0 + c0, &carry, &n0);

   // Compute n1 and whatever is carried.
   Split (carry + w1 + x0 + y0 + c1, &carry, &n1);

   // Compute c0 and whatever is carried.
   Split (carry + x1 + y1 + z0, &carry, &c0);

   // Compute c1; nothing is carried.
   c1 = z1 + carry;

   // Re-assemble n1 and n0.
   N = (n1 << 16) + n0;

   return ((N + 0.5) / 4294967296.0);

}

// This function computes 16-bit numbers x1 and x0 such that x = x1 * 2^16 + x0.
// Pointers *x1 and *x0 are used because a function can return only one value.
void Split (unsigned int x, unsigned int *x1, unsigned int *x0) {

   // First 16 bits of x.
   *x0 = x & 0xffff;

   // Second 16 bits of x.
   *x1 = x >> 16;

   return;

}

// This is the Mersenne Twister Temper function.
unsigned int Temper (unsigned int N) {

   N ^= (N >> 11);
   N ^= (N << 7) & 0x9d2c5680;
   N ^= (N << 15) & 0xefc60000;
   N ^= (N >> 18);

   return N;

}











