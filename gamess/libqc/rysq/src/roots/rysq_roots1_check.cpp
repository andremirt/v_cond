#include <math.h>

//#include "rysq_roots1.h"
#include "rysq_roots1.h.allada"

#include <stdlib.h>
#include <time.h>
#include <limits>
#include <iostream>


void Rysq_roots1_check() {
    
    double X_range[8] = {33.0, 15.0, 10.0, 5.0, 3.0, 1.0, 3e-7, 0.0};
    
    srand((unsigned)time(NULL));

    for(int i = 0; i < 14; ++i) {
      double t2[1] = {0.0};
      double t2_chk[1] = {0.0};
      double W[1] = {0.0};
      double W_chk[1]= {0.0};
      
      double interval = (i/2 > 0) ? X_range[i/2-1] - X_range[i/2] :  1000.0;
      double X = X_range[i/2] + interval*((double)rand()/RAND_MAX);
      
      
      // 	if (X > 33.) {
      // 	    rysq_roots1_x33toINF(X, t2, W);
      // 	    Rysq_roots1_x33toINF(X, t2_chk, W_chk);
      // 	}
      // 	else if(X > 15.) {
      // 	    rysq_roots1_x15to33(X, t2, W);
      // 	    Rysq_roots1_x15to33(X, t2_chk, W_chk);
      // 	}
      // 	else if (X > 10.) {
      // 	    rysq_roots1_x10to15(X, t2, W);
      // 	    Rysq_roots1_x10to15(X, t2_chk, W_chk);
      // 	}
      //	else if (X > 5.) {
      // 	    rysq_roots1_x5to10(X, t2, W);
      // 	    Rysq_roots1_x5to10(X, t2_chk, W_chk);
      // 	}
      // 	else if (X > 3.) {
      // 	    rysq_roots1_x3to5(X, t2, W);
      // 	    Rysq_roots1_x3to5(X, t2_chk, W_chk);
      // 	}
      // 	else if (X > 1.) {
      // 	    rysq_roots1_x1to3(X, t2, W);	 
      // 	    Rysq_roots1_x1to3(X, t2_chk, W_chk);
      // 	}
      // 	else if (X > 3e-7) {
      // 	    rysq_roots1_x0to1(X, t2, W);
      // 	    Rysq_roots1_x0to1(X, t2_chk, W_chk);
      // 	}
      // 	/*     X IS APPROXIMATELY ZERO.         NROOTS=1 */
      // 	else {
      // 	    t2[0] = .5 - X / 5.;
      // 	    W[0] = 1. - X / 3.;
      // 	}
      
      Rysq_roots1(X, t2, W);
      Rysq_roots1(X, t2_chk, W_chk);
      
      printf("X = %f, t2 = %f, t2_check = %f, error = %f ", X, t2[0], t2_chk[0], (t2[0]-t2_chk[0])/t2[0]);
      printf("W = %f, W_check = %f, error = %f\n",W[0], W_chk[0], (W[0]-W_chk[0])/W[0]);
      
      // 	for(int i = 0; i < 1; ++i) {
      // 	    t2[i] = t2[i]/(1.0 + t2[i]);
      // 	}
    }
    return;
}

int main() {
  
    Rysq_roots1_check();
    return 0;
}
