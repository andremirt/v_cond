#include <math.h>

//#include "rysq_roots1.h"
#include "rysq_roots3.h.allada"
#include <stdlib.h>
#include <time.h>
#include <limits>
#include <iostream>

void Rysq_roots3_check() {
    
    double X_range[9] = {47.0, 33.0, 15.0, 10.0, 5.0, 3.0, 1.0, 3e-7, 0.0};
    srand((unsigned)time(NULL));
    
    for(int i = 0; i < 16; ++i) {
	double t2[3] = {0.0, 0.0, 0.0};
	double t2_chk[3] = {0.0, 0.0, 0.0};
	double W[3] = {0.0, 0.0, 0.0};
	double W_chk[3]= {0.0, 0.0, 0.0};
	
	double interval = (i/2 > 0) ? X_range[i/2-1] - X_range[i/2] :  1000.0;
	double X = X_range[i/2] + interval*((double)rand()/RAND_MAX);
	
	if (X > 47.) {
	    rysq_roots3_x47toINF(X, t2, W);
	}
	else if (X > 33.) {
	    rysq_roots3_x33to47(X, t2, W);
	}
	else if (X > 20.) {
	    rysq_roots3_x20to33(X, t2, W);
	}
	else if (X > 15.) {
	    Ghondo_roots3_x15to20(X, t2, W);
	    Rysq_roots3_x15to20(X, t2_chk, W_chk);
	}
	else if (X > 10.) {
	    Ghondo_roots3_x10to15(X, t2, W);
	    Rysq_roots3_x10to15(X, t2_chk, W_chk);
		}
	else if (X > 5.) {
	    Ghondo_roots3_x5to10(X, t2, W);
	    Rysq_roots3_x5to10(X, t2_chk, W_chk);
		}
	else if (X > 3.) {
	    Ghondo_roots3_x3to5(X, t2, W);
	    Rysq_roots3_x3to5(X, t2_chk, W_chk);
	}
	else if (X > 1.) {
	    Ghondo_roots3_x1to3(X, t2, W);
	    Rysq_roots3_x1to3(X, t2_chk, W_chk);
		}
	else if (X > 3e-7) {
	    Ghondo_roots3_x0to1(X, t2, W);
	    Rysq_roots3_x0to1(X, t2_chk, W_chk);
	}
	else {
	    /*     X IS APPROXIMATELY ZERO.         NROOTS=3 */
	    t2[0] = .0603769246832797 - X * .00928875764357368;
	    t2[1] = .776823355931043 - X * .119511285527878;
	    t2[2] = 6.66279971938567 - X * 1.02504611068957;
	    W[0] = .467913934572691 - X * .0564876917232519;
	    W[1] = .360761573048137 - X * .149077186455208;
	    W[2] = .171324492379169 - X * .127768455150979;
	}
	
	printf ("X = %f\n", X);
	printf("t2[0] = %f, t2_check[0] = %f, error = %f\n", t2[0], t2_chk[0], (t2[0]-t2_chk[0])/t2[0]);
	printf("t2[1] = %f, t2_check[1] = %f, error = %f\n ", t2[1], t2_chk[1], (t2[1]-t2_chk[1])/t2[1]);
	printf("t2[2] = %f, t2_check[2] = %f, error = %f\n ", t2[2], t2_chk[2], (t2[2]-t2_chk[2])/t2[2]);
	printf("W[0] = %f, W_check[0] = %f, error = %f\n", W[0], W_chk[0], (W[0]-W_chk[0])/W[0]);
	printf("W[1] = %f, W_check[1] = %f, error = %f\n", W[1], W_chk[1], (W[1]-W_chk[1])/W[1]);
	printf("W[2] = %f, W_check[2] = %f, error = %f\n", W[2], W_chk[2], (W[2]-W_chk[2])/W[2]);
	
	//       Rysq_roots1(X, t2, W);
	//       Rysq_roots1(X, t2_chk, W_chk);
    }
    return;
}

int main() {
    Rysq_roots3_check();
    return 0;
}
