/**
   @file
   @brief Rys quadrature roots and weights implementation
*/

#include "rysq_roots_flymake.h"
#include "opq.h"
#include "rysq_roots0_flymake.h"
#include "rysq_roots1_flymake.h"
#include "rysq_roots2.h"
#include "rysq_roots3.h"
#include "rysq_roots4.h"
#include "rysq_roots5.h"
#include <stdlib.h>
#include <assert.h>

static bool initialized = false;



int Rysq_roots_init() {
    rysq::stieltjes_init();
    initialized = true;
    return 0;
}

int Rysq_roots_finalize() {
    rysq::stieltjes_finalize();
    initialized = false;
    return 0;
}

void Rysq_roots(int N, double X, double *t2, double *W) {
    assert(0<=N && N<=13);

    if (N == 0) {
	W[0] = Rysq_roots0(X);
    }
    else if (N == 1) {
	Rysq_roots1(X, t2, W);
    }
    else if (N == 2) {
	Rysq_roots2(X, t2, W);
    }
    else if (N == 3) rysq::roots<3>(X, t2, W);
    else if (N == 4) rysq::roots<4>(X, t2, W);
    else if (N == 5) rysq::roots<5>(X, t2, W);
    else {
	assert(initialized);
	//rysq::stieltjes(N, X, t2, W);
    }

    return;
}

double *rysq::STIELTJES_X[8]; 

double *rysq::STIELTJES_W[8]; 

/**
   @brief Initializes auxiliary Stieltjes quadratures
*/
void rysq::stieltjes_init() {
    double a[55], b[55];

    for(int i = 0; i < 8; ++i) {
	int N = STIELTJES_N[i]; 

	STIELTJES_X[i] = (double*)malloc(N*sizeof(double));
	STIELTJES_W[i] = (double*)malloc(N*sizeof(double));

	a[0] = 0.5;
	b[0] = 1.0; 

	for(int j = 1; j < N; ++j) {
	    a[j] = 0.5;
	    b[j] = 0.25/(4.0 - (1.0/(j*j)));
	}

	OPQ_coeff(N, a, b, STIELTJES_X[i], STIELTJES_W[i]); 
    }
}

/**
   @brief Finalizes auliliary Stieltjes quadratures
*/
void rysq::stieltjes_finalize() {
    
    for(int i = 0; i < 8; ++i) {

	if(STIELTJES_X[i]) {
	    free(STIELTJES_X[i]);
	    STIELTJES_X[i] = NULL;
	}

	if(STIELTJES_W[i]) {
	    free(STIELTJES_W[i]);
	    STIELTJES_W[i] = NULL;
	}

    }
}

