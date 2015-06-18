#include <math.h>

//#include "rysq_roots1.h"
#include "rysq_roots2.h.allada"
#include <stdlib.h>
#include <time.h>
#include <limits>
#include <iostream>

void Rysq_roots2_check() {
    
    double X_range[9] = {40.0, 33.0, 15.0, 10.0, 5.0, 3.0, 1.0, 3e-7, 0.0};
    srand((unsigned)time(NULL));
    
    for(int i = 0; i < 16; ++i) {
	double t2[2] = {0.0, 0.0};
	double t2_chk[2] = {0.0, 0.0};
	double W[2] = {0.0, 0.0};
	double W_chk[2]= {0.0, 0.0};
	
	double interval = (i/2 > 0) ? X_range[i/2-1] - X_range[i/2] :  1000.0;
	double X = X_range[i/2] + interval*((double)rand()/RAND_MAX);
	
	if (X > 40.) {
	    Ghondo_roots2_x40toINF(X, t2, W);
	}
	else if (X > 33.) {
	    Ghondo_roots2_x33to40(X, t2, W);
	    Rysq_roots2_x33to40(X, t2_chk, W_chk);
	}
	else if (X > 15.) {
	    Ghondo_roots2_x15to33(X, t2, W);
	    Rysq_roots2_x15to33(X, t2_chk, W_chk);
	}  
	else if (X > 10.) {
	    Ghondo_roots2_x10to15(X, t2, W);
	    Rysq_roots2_x10to15(X, t2_chk, W_chk);
	}
	else if (X > 5.) {
	    Ghondo_roots2_x5to10(X, t2, W);
	    Rysq_roots2_x5to10(X, t2_chk, W_chk);
	}
	else if (X > 3.) {
	    Ghondo_roots2_x3to5(X, t2, W);
	    Rysq_roots2_x3to5(X, t2_chk, W_chk);
	}
	else if (X > 1.) {
	    Ghondo_roots2_x1to3(X, t2, W);
	    Rysq_roots2_x1to3(X, t2_chk, W_chk);
	}
	else if (X > 3.0e-7) {
	    Ghondo_roots2_x0to1(X, t2, W);
	    Rysq_roots2_x0to1(X, t2_chk, W_chk);
	}
	else {
	    /*     X IS APPROXIMATELY ZERO.         NROOTS=2 */
	    t2[0] = .130693606237085 - X * .0290430236082028;
	    t2[1] = 2.86930639376291 - X * .637623643058102;
	    W[0] = .652145154862545 - X * .122713621927067;
	    W[1] = .347854845137453 - X * .210619711404725;
	}
	
    	for(int i = 0; i < 2; ++i) {
    	    t2[i] = t2[i]/(1.0 + t2[i]);
    	    t2_chk[i] = t2_chk[i]/(1.0 + t2_chk[i]);
    	}
	
	printf ("X = %f\n", X);
	printf("t2[0] = %f, t2_check[0] = %f, error = %f\n", t2[0], t2_chk[0], (t2[0]-t2_chk[0])/t2[0]);
	printf("t2[1] = %f, t2_check[1] = %f, error = %f\n ", t2[1], t2_chk[1], (t2[1]-t2_chk[1])/t2[1]);
	printf("W[0] = %f, W_check[0] = %f, error = %f\n", W[0], W_chk[0], (W[0]-W_chk[0])/W[0]);
	printf("W[1] = %f, W_check[1] = %f, error = %f\n", W[1], W_chk[1], (W[1]-W_chk[1])/W[1]);
	
    }
    
    
    //       Rysq_roots1(X, t2, W);
    //       Rysq_roots1(X, t2_chk, W_chk);
    
    return;
}

int main() {
  
    Rysq_roots2_check();
    return 0;
}
