#include <iostream>
#include "Functions.h"

int main () {
    
   printf ("The size of an unsigned...\n");
   printf (" - short int is %d bytes\n", sizeof (unsigned short int));
   printf (" - int is %d bytes\n", sizeof (unsigned int));
   printf (" - long int is %d bytes\n", sizeof (unsigned long int));

   int n;
   Time ();
   for (n = 1; n <= 500000000; n++) MTUniform(0);

   printf ("Computations took %5.2f seconds.\n", Time());

   Pause ();

}   







