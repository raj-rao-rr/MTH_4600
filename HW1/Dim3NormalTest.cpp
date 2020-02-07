
#include "Functions.h"

// MTUniform is the Mersenne twister, MWCUniform is the multiply with carry
//    generator that we discussed in class.


double MWCUniform (unsigned int);
unsigned int Temper (unsigned int);
int power (int, int);
int Uniform(int);
int test(int arr[]);


int power(int x, int a) {
    int val = 1;

    // iterates to the power provided in order to determine the expression x^a
    for (int i = 0; i < a; ++i) {
        val *= x;

    return val;
}
    

int Uniform(int w) {
    // determines which uniform random number generator is to be used either MTUniform (Mersenne twister) or MWCUniform (multiply with carry)
    if (w == 1) {
        return MTUniform(0);
    }
    else if (w == 2) {
        return MWCUniform(0);
    }
}


int test(int U[]) {
    double delta, M = 0;
    int pw = 3, c=0;

    // split the unit interval into 40 disjoint subintervals of equal length 
    delta = 1 / 40;

    // k iterates from 0 < k < 40 for each disjoint interval
    for (int k = 0; k < 40; ++k) {
        // creating the bounds for the restricted interval I, (k*delta, (k+1)*delta)
        lower = k * delta;
        upper = (k + 1) * delta;

        // if we have that our uniform random variable U[c] at index (c) is in the restricted interval I denoted by the lower/upper levels  
        if ((U[c] > lower) && (U[c] <= upper) && pw > 0) {
            // if the condition is met we form our M through the sumproduct of k's and 40 raised to the powers [2,1,0] in that order
            M += k * power(40, pw - 1);
            --pw;
            ++c;
        }
    }
    return M;
}


int main() {
    double m, w, n, U[3];
    int* X;

    // Allocate space for X[1],... X[64000] and initialize each to 0.
    X = (int*)calloc(64000, sizeof(int));

    // allow the user to select their specific number generator 
    std::cout << "Which random number generator would you like to implement";
    std::cout << "(1) MTUniform";
    std::cout << "(2) MWCUniform";
    std::cin >> w;

    std::cout << "How many simulations would you like to run (suggested 10 million)";
    std::cin >> n;
    

    // running n simulations (suggested 10, 50, 100 million) to test normalacy
    for (int i = 0; i < n; ++i) {

        for (int j = 0; j < 3; ++j) {
            U[j] = Uniform(w);
        }

        // Compute the value of M provided the testing instructions
        m = test(U);

        // Increment the appropriate X by 1.
        X[m] ++;
    
    }

    // Below code adopted from C. Douglas Howard.

    // Now each X[m] should be Binomial (n, p), where p = 1/40320.
    p = 1.0 / 64000.0;
    q = 1.0 - p;
    mu = n * p;
    sigma = sqrt(n * p * q);

    // If n is large, (X[m] - mu)/sigma should be approximately Normal(0,1).
    // Normalize and add to a normal histogram.
    for (m = 1; m <= 64000; ++m) {
        Z = (X[m] - mu) / sigma;
        NormalHistogram(Z, 40, 0);
    }

    // Create the TeX files for viewing.
    NormalHistogram(0, 40, 1);
    
    std::cout << "The simulation is complete...";
    
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











