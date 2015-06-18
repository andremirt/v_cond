#include <math.h>

//#include "rysq_roots0.h"       	
//#include "rysq_roots0.h.mine"
#include "rysq_roots0.h.allada"

#include <stdlib.h>
#include <time.h>
#include <limits>
#include <iostream>


void Rysq_roots0_check() {
  
  double X_range[8] = {33.0, 15.0, 10.0, 5.0, 3.0, 1.0, 3e-7, 0.0};
  double t2 = 0.0; 
  double t2_chk = 0.0;
  
  srand((unsigned)time(NULL));
  
  for(int i = 0; i < 14; ++i) {
    double interval = (i/2 > 0) ? X_range[i/2-1] - X_range[i/2] :  1000.0;
    double X = X_range[i/2] + interval*((double)rand()/RAND_MAX);
    
    t2 = Ghondo_Rysq_roots0(X);
    t2_chk = Rysq_roots0(X);

//     if (X > 33.) {
//       t2 = Ghondo_roots0_x33toInf(X);
//       t2_chk = Rysq_roots0_x33toInf(X);
//     }
//     else if(X > 15.) {
//       t2 = Ghondo_roots0_x15to33(X);
//       t2_chk = Rysq_roots0_x15to33(X);
//     }
//     else if (X > 10.) {
//       t2 = Ghondo_roots0_x10to15(X);
//       t2_chk = Rysq_roots0_x10to15(X);
//     }
//     else if (X > 5.) {
//       t2 = Ghondo_roots0_x5to10(X);
//       t2_chk = Rysq_roots0_x5to10(X);
//     }
//     else if (X > 3.) {
//       t2 = Ghondo_roots0_x3to5(X);
//       t2_chk = Rysq_roots0_x3to5(X);
      
//     }
//     else if (X > 1.) {
//       t2 = Ghondo_roots0_x1to3(X);
//       t2_chk = Rysq_roots0_x1to3(X);
//     }
//     else if (X > 3e-7) {
//       t2 = Ghondo_roots0_x0to1(X);
//       t2_chk = Rysq_roots0_x0to1(X);
//     }
//     /*     X IS APPROXIMATELY ZERO.         NROOTS=1 */
//     else {
//       t2 = (1. - X / 3.);
//     }

    printf("X = %f, t2 = %f, t2_chk = %f, error = %f \n", X, t2, t2_chk, (t2 - t2_chk)/t2);
    
  }
  return;
}

int main() {

    Rysq_roots0_check();
    return 0;
}
