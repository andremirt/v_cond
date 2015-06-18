/**
   @file
*/


#include "eri/rysq_eri_int2d.h"
#include "cuda/curysq_const.h"
#include "cuda/curysq_util.h"
#include "vec.h"
#include "util.h"
#include <math.h>
#include "roots/rysq_roots2.h"
#include <stdio.h>
#include <iostream>
#include <algorithm>

#define pow2(x) ((x)*(x))
#define pow3(x) ((x)*(x)*(x))

extern __shared__ double shmem[];

typedef struct {
    double rbra[6];
    double rket[6];
    double rij[3];
    double rkl[3];
    double rij2;
    double rkl2;
    double *eij, *ekl;
    double *rA, *rB;
    double *A, *A1, *B, *B1;
    double *rAi, *rBk;
} sh2_t;

/**
   @brief Set up shared memory pointers
*/
__device__ void sharedInit(bool braRecur, bool ketRecur, int Kbra, int Kket, 
			   double2 *dev_braket, double *dev_rbra, double *dev_rket,
			   sh2_t &sh) {

    short rank = ctaRank();
    short bsize = ctaSize();

    // load bra and ket centers
    for(int i = rank; i < 6; i += bsize) {
	sh.rbra[i] = dev_rbra[i];
	sh.rket[i] = dev_rket[i];
    }
    __syncthreads();    
   
    // compute bra and ket distances
    for(int i = rank; i < 3; i += bsize) {
	sh.rij[i] = sh.rbra[i] - sh.rbra[i+3];
	sh.rkl[i] = sh.rket[i] - sh.rket[i+3];
    }
    __syncthreads();


    if(rank == 0) {
	sh.rij2 = 0.0;
	sh.rkl2 = 0.0;
	for(int i = 0; i < 3; ++i) {
	    sh.rij2 += pow2(sh.rij[i]);
	    sh.rkl2 += pow2(sh.rkl[i]);
	}
	double *sh_next = shmem;
	
	// computed bra values
	sh.rA = sh_next;  sh_next += 3*Kbra;
	sh.eij = sh_next; sh_next += Kbra;

	// computed ket values
	sh.rB = sh_next;  sh_next += 3*Kket;
	sh.ekl = sh_next; sh_next += Kket;

	// recurrence variables
	sh.A = sh_next;  sh_next += Kbra;
	sh.A1 = sh_next; sh_next += Kbra;
	sh.B = sh_next;  sh_next += Kket;
	sh.B1 = sh_next; sh_next += Kket;
	
	if(braRecur)
	    sh.rAi = sh_next; sh_next += 3*Kbra;
	
	if(ketRecur)
	    sh.rBk = sh_next; sh_next += 3*Kket;
    }
    __syncthreads();

    // compute bra values
    for(int Kij = rank; Kij < Kbra; Kij += bsize) {
	double ai = dev_braket[Kij].x;
	double aj = dev_braket[Kij].y;
	double A = ai + aj;
	double A1 = 1.0/A;

	sh.eij[Kij] = exp(-ai*aj*A1*sh.rij2);
	for(int i=0; i<3; ++i) {
	    double rA = A1*(ai*sh.rbra[i] + aj*sh.rbra[i+3]);
	    if(braRecur)
		sh.rAi[i+Kij*3] = rA - sh.rbra[i];
	    sh.rA[i+Kij*3] = rA;
	}
	sh.A[Kij] = A;
        sh.A1[Kij] = A1;
    }

    // compute ket values
    for(int Kkl = rank; Kkl < Kket; Kkl += bsize) {
	double ak = dev_braket[Kkl+Kbra].x;
	double al = dev_braket[Kkl+Kbra].y;
	double B = ak + al;
	double B1 = 1.0/B;

	sh.ekl[Kkl] = exp(-ak*al*B1*sh.rkl2);
	for(int i = 0; i < 3; ++i) {
	    double rB = B1*(ak*sh.rket[i] + al*sh.rket[i+3]);
	    if(ketRecur)
		sh.rBk[i+Kkl*3] = rB - sh.rket[i];
	    sh.rB[i+Kkl*3] = rB;
	}
	sh.B[Kkl] = B;
	sh.B1[Kkl] = B1;
    }
    __syncthreads();
}



/**
   @brief
   @param a
   @param b
*/
void cuRysq_packExp_r2(Rysq_shell_t a, Rysq_shell_t b,
		       Rysq_shell_t c, Rysq_shell_t d,
		       double2 *exp) {
    
    int Kbra = a.K*b.K;

    for(int j = 0, ij = 0; j < b.K; ++j) {
	for(int i = 0; i < a.K; ++i, ++ij) {
	    exp[ij].x = a.a[i];
	    exp[ij].y = b.a[j];
	}
    }

    for(int l = 0, kl = 0; l < d.K; ++l) {
	for(int k = 0; k < c.K; ++k, ++kl) {
	    exp[kl+Kbra].x = c.a[k];
	    exp[kl+Kbra].y = d.a[l];
	}
    }
}


/**
   @brief
   @param flags
*/
__global__ void cuRysq_eri_dsss(int flags, double tole, 
				double2 *dev_AB2rho,
				int Kbra, int Kket, double2 *dev_braket,
				double *dev_rbra, double *dev_rket,
				double scale, double *dev_I) {
    
    short rank = ctaRank();
    short bsize = ctaSize();

    __shared__ sh2_t sh;
    sharedInit(true, false, Kbra, Kket, dev_braket, dev_rbra, dev_rket, sh);

    double q[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    // compute contractions
    for(unsigned short int Kkl = threadIdx.z; Kkl < Kket; Kkl += blockDim.z) {

	//the bra Kij contractions are mapped to threadIdx.y
	const unsigned short Kij = threadIdx.y;
	const unsigned short K = Kij + Kkl*Kbra;

	double AB2 = dev_AB2rho[K].x;
	double rho = dev_AB2rho[K].y;
	double CAB2e = AB2*sh.eij[Kij]*sh.ekl[Kkl];
	if (fabs(CAB2e) < tole) continue;
	
	double X = 0.0;
	for(int i = 0; i < 3; ++i) {
	    double rAB = (sh.rA[i+Kij*3]-sh.rB[i+Kkl*3]);
	    X += pow2(rAB);
	}		
	X *= (sh.A[Kij]*sh.B[Kkl])*rho;

	double t2[2], W[2];
	Rysq_roots2(X, t2, W);

	t2[0] *= rho;
	t2[1] *= rho;

	double WBt[3];
	WBt[0] = W[0] + W[1];
	WBt[1] = sh.B[Kkl]*(W[0]*t2[0] + W[1]*t2[1]);
	WBt[2] = pow2(sh.B[Kkl])*(W[0]*t2[0]*t2[0] + W[1]*t2[1]*t2[1]);
	
	double rAB = (sh.rA[Kij*3] - sh.rB[Kkl*3]);
	double rAi = sh.rAi[Kij*3];
	double t1 = 0.5*sh.A1[Kij]*(WBt[0] - WBt[1]);	

	//last integrals have to be swapped
	for(int i = 0; i < 3; ++i) {
	    double x1 = WBt[2]*rAB - WBt[1]*rAi;
	    double x2 = WBt[0]*rAi - WBt[1]*rAB;

	    q[i] += CAB2e*(rAB*x1 + rAi*x2 + t1);

	    int j = (i+1)%3;
	    rAB = (sh.rA[j+Kij*3] - sh.rB[j+Kkl*3]);
	    rAi = sh.rAi[j+Kij*3];

	    q[i+3] += CAB2e*(rAB*x1 + rAi*x2);
 	}
    }
    __syncthreads();

    // reduce values
    q[0] = ctaReduce(q[0],  shmem);
    q[1] = ctaReduce(q[1],  shmem);
    q[2] = ctaReduce(q[2],  shmem);
    q[3] = ctaReduce(q[3],  shmem);
    q[4] = ctaReduce(q[4],  shmem);
    q[5] = ctaReduce(q[5],  shmem);

    // put values into shared memory
    if(rank == 0) {
	for(int i = 0; i < 6; ++i) {
	    shmem[i] = q[i];
	}
	// last two integrals have to be swapped
	swap(shmem[4], shmem[5]);
    }
    __syncthreads();

    // scale and put values into global memory
    for(int i = rank; i < 6; i += bsize) {
	dev_I[i] = scale*SQRT_4PI5*CURYSQ_NORMAL[i+4]*shmem[i];
    }
}

/**
   @brief
   @param flags
*/
__global__ void cuRysq_eri_fsss(int flags, double tole, 
				double2 *dev_AB2rho,
				int Kbra, int Kket, double2 *dev_braket,
				double *dev_rbra, double *dev_rket,
				double scale, double *dev_I) {
    
    short rank = ctaRank();
    short bsize = ctaSize();

    __shared__ sh2_t sh;
    sharedInit(true, false, Kbra, Kket, dev_braket, dev_rbra, dev_rket, sh);

    double q[10];
    for(int r = 0; r < 10; ++r) {
	q[r] = 0.0;
    }
    
    // compute contractions
    for(unsigned short int Kkl = threadIdx.z; Kkl < Kket; Kkl += blockDim.z) {
	//the bra Kij contractions are mapped to threadIdx.y
	const unsigned short Kij = threadIdx.y;
	const unsigned short K = Kij + Kkl*Kbra;

	double AB2 = dev_AB2rho[K].x;
	double rho = dev_AB2rho[K].y;
	double CAB2e = AB2*sh.eij[Kij]*sh.ekl[Kkl];
	if (fabs(CAB2e) < tole) continue;
	
	double X = 0.0;
	for(int i = 0; i < 3; ++i) {
	    double rAB = (sh.rA[i+Kij*3]-sh.rB[i+Kkl*3]);
	    X += pow2(rAB);
	}		
	X *= (sh.A[Kij]*sh.B[Kkl])*rho;
	double t2[2], W[2];
	Rysq_roots2(X, t2, W);

	t2[0] *= rho;
	t2[1] *= rho;

	// create root values into memory
	double WBt[4];
	WBt[0] = W[0] + W[1];
	WBt[1] = sh.B[Kkl]*(W[0]*t2[0] + W[1]*t2[1]);
	WBt[2] = pow2(sh.B[Kkl])*(W[0]*t2[0]*t2[0] + W[1]*t2[1]*t2[1]);
	WBt[3] = pow3(sh.B[Kkl])*(W[0]*pow3(t2[0]) + W[1]*pow3(t2[1]));
	
	// solve recurrence
	//(3,0,0)
	//(0,3,0)
	//(0,0,3)
	//(2,1,0)
	//(0,2,1)
	//(1,0,2)
	//(2,0,1)
	//(1,2,0)
	//(0,1,2)
	double rAB = (sh.rA[Kij*3]-sh.rB[Kkl*3]);
	double rAi = sh.rAi[Kij*3];	
	for(int i = 0; i < 3; ++i) {
	    double x1 = 0.5*sh.A1[Kij]*(WBt[0] - WBt[1]);
	    double x2 = WBt[0]*pow2(rAi) - WBt[1]*2*rAB*rAi + WBt[2]*pow2(rAB);
	    double x3 = (WBt[2] - WBt[1])*0.5*sh.A1[Kij];
	    double x4 = WBt[2]*2*rAB*rAi - WBt[1]*pow2(rAi) - WBt[3]*pow2(rAB);
	    
	    q[i] += CAB2e*(x1*(3*rAi) + x2*(rAi) + x3*(3*rAB) + x4*(rAB));
	    
	    int j = (i+2)%3;
	    rAB = (sh.rA[j+Kij*3]-sh.rB[j+Kkl*3]);
	    rAi = sh.rAi[j+Kij*3];
	    q[i+6] += CAB2e*(x1*(rAi) + x2*(rAi) + x3*(rAB) + x4*(rAB));
	    
	    j = (i+1)%3;
	    rAB = (sh.rA[j+Kij*3]-sh.rB[j+Kkl*3]);
	    rAi = sh.rAi[j+Kij*3];
	    q[i+3] += CAB2e*(x1*(rAi) + x2*(rAi) + x3*(rAB) + x4*(rAB));
	}
	//rAB is xAB
	double yAB = (sh.rA[1+Kij*3]-sh.rB[1+Kkl*3]);
	double zAB = (sh.rA[2+Kij*3]-sh.rB[2+Kkl*3]);
	//(1,1,1)
	q[9] += CAB2e*(WBt[0]*(rAi*sh.rAi[1+Kij*3]*sh.rAi[2+Kij*3]) -
		       WBt[1]*(rAB*sh.rAi[1+Kij*3]*sh.rAi[2+Kij*3] + rAi*yAB*sh.rAi[2+Kij*3] + rAi*sh.rAi[1+Kij*3]*zAB) +
		       WBt[2]*(rAi*yAB*zAB + rAB*sh.rAi[1+Kij*3]*zAB + rAB*yAB*sh.rAi[2+Kij*3]) -
		       WBt[3]*(rAB*yAB*zAB));

    }
    __syncthreads();
	
    // reduce values
    q[0] = ctaReduce(q[0],  shmem);
    q[1] = ctaReduce(q[1],  shmem);
    q[2] = ctaReduce(q[2],  shmem);
    q[3] = ctaReduce(q[3],  shmem);
    q[4] = ctaReduce(q[4],  shmem);
    q[5] = ctaReduce(q[5],  shmem);
    q[6] = ctaReduce(q[6],  shmem);
    q[7] = ctaReduce(q[7],  shmem);
    q[8] = ctaReduce(q[8],  shmem);
    q[9] = ctaReduce(q[9],  shmem);

    // put values into shared memory
    if(rank == 0) {
	for(int i = 0; i < 10; ++i) {
	    shmem[i] = q[i];
	}
	//put integrals into correct order
	//(3,0,0)
	//(0,3,0)
	//(0,0,3)
	//(2,1,0)
	//(2,0,1)
	//(1,2,0)
	//(1,0,2)
	//(0,2,1)
	//(0,1,2)
	swap(shmem[4], shmem[6]);
	swap(shmem[6], shmem[7]);
	swap(shmem[5], shmem[6]);
    }
    __syncthreads();

    // scale and put values into global memory
    for(int i = rank; i < 10; i += bsize) {
	dev_I[i] = scale*SQRT_4PI5*shmem[i];
    }
}


/**
   @brief
   @param flags
*/
__global__ void cuRysq_eri_psps(int flags, double tole, 
				double2 *dev_AB2rho,
				int Kbra, int Kket, double2 *dev_braket,
				double *dev_rbra, double *dev_rket,
				double scale, double *dev_I) {
    
    short rank = ctaRank();
    short bsize = ctaSize();

    __shared__ sh2_t sh;
    sharedInit(true, true, Kbra, Kket, dev_braket, dev_rbra, dev_rket, sh);

    // initialize q values
    double q[9];
    for(int r = 0; r < 9; ++r) {
	q[r] = 0;
    }
    
    // compute contractions
    for(unsigned short int Kkl = threadIdx.z; Kkl < Kket; Kkl += blockDim.z) {
	//the bra Kij contractions are mapped to threadIdx.y
	const unsigned short Kij = threadIdx.y;
	const unsigned short K = Kij + Kkl*Kbra;

	double AB2 = dev_AB2rho[K].x;
	double rho = dev_AB2rho[K].y;
	double CAB2e = AB2*sh.eij[Kij]*sh.ekl[Kkl];
	//if (fabs(CAB2e) < tole) continue;
	
	double X = 0.0;
	for(int i = 0; i < 3; ++i) {
	    double rAB = (sh.rA[i+Kij*3]-sh.rB[i+Kkl*3]);
	    X += pow2(rAB);
	}		
	X *= (sh.A[Kij]*sh.B[Kkl])*rho;
	double t2[2], W[2];
	Rysq_roots2(X, t2, W);

	t2[0] *= rho;
	t2[1] *= rho;

	// create root values into shared memory
	double Wt[3];
	Wt[0] = (W[0] + W[1]);
	Wt[1] = (W[0]*t2[0] + W[1]*t2[1]);
	Wt[2] = (W[0]*t2[0]*t2[0] + W[1]*t2[1]*t2[1]);
	
	// solve recurrence
	for(int i = 0; i < 3; ++i) {
	    double rAB = (sh.rA[i+Kij*3]-sh.rB[i+Kkl*3]);
	    q[i] += CAB2e*(Wt[0]*(sh.rAi[i+Kij*3]*sh.rBk[i+Kkl*3]) +
			   Wt[1]*(0.5 + sh.A[Kij]*rAB*sh.rAi[i+Kij*3] - sh.B[Kkl]*rAB*sh.rBk[i+Kkl*3]) -
			   Wt[2]*(sh.A[Kij]*sh.B[Kkl]*pow2(rAB)));
	}

	
	for(int bk = 0, r = 3; bk < 3; ++bk) {
	    double rABbk = (sh.rA[bk+Kij*3]-sh.rB[bk+Kkl*3]);
	    for(int ai = 0; ai < 3; ++ai) {
		if(ai == bk) continue;
		double rABai = (sh.rA[ai+Kij*3]-sh.rB[ai+Kkl*3]);
		q[r] += CAB2e*(Wt[0]*(sh.rAi[ai+Kij*3]*sh.rBk[bk+3*Kkl]) +
			       Wt[1]*(sh.A[Kij]*sh.rAi[ai+Kij*3]*rABbk - sh.B[Kkl]*rABai*sh.rBk[bk+Kkl*3]) -
			       Wt[2]*(sh.A[Kij]*sh.B[Kkl]*rABai*rABbk));
		++r;
	    }
	}
    }

    __syncthreads();

    // reduction summation
    __shared__ double I[9];

    //TODO: reduction needs to be in the following forms
    //CTA reduce needs to be implemented
    I [0] =ctaReduce (q [0],shmem);
    __syncthreads();
    // ...;
    __syncthreads();
    I [8] =ctaReduce (q [8],shmem);
    __syncthreads();

    for(int i = rank; i < 9; i += bsize) {
	dev_I[i] = scale*SQRT_4PI5*I[i];
    }
}


/**
   @brief
*/
int cuRysq_eri1_r2(int flags, double tol,
		   Rysq_shell_t a, double *ri,
		   Rysq_shell_t b, double *rj,
		   Rysq_shell_t c, double *rk,
		   Rysq_shell_t d, double *rl,
		   double scale, double *I) {

    //int L = a.L + b.L + c.L + d.L;
    int mask = mask(0,1,2,3);
    //     if(L > 1) {
    // 	if((a.L + b.L + c.L ) == 0) {
    // 	    mask = mask(3, 2, 1, 0);
    // 	} else if((a.L + c.L) == 0) {
    // 	    mask = mask(1, 0, 3, 2);
    //  	} else if((a.L + b.L) == 0) {
    //  	    mask = mask(2, 3, 0, 1);
    // 	} else if(c.L == 0) {
    // 	    mask = mask(0, 1, 3, 2);
    // 	} else if(a.L == 0) {
    // 	    mask = mask(1, 0, 2, 3);
    // 	}
    //     }

    double *r1 = ri; 
    double *r2 = rj;
    double *r3 = rk; 
    double *r4 = rl;
    shuffle(a,b,c,d,mask);
    shuffle(r1,r2,r3,r4,mask);

    bool sp = false;

    int Kbra = a.K*b.K;
    int Kket = c.K*d.K;

    double2 braket[Kbra+Kket];
    cuRysq_packExp_r2(a, b, c, d, braket);

    //AB2 includes contrations
    double2 AB2rho[Kbra*Kket];
    //Csp is used to create s contraction coefficients
    double Csp[Kket];

    for(int l = 0, ijkl = 0, kl = 0; l < d.K; ++l) {
	for(int k = 0; k < c.K; ++k, ++kl) {
	    double B = braket[kl+Kbra].x+braket[kl+Kbra].y;
	    double ckl;
	    if(sp) {
		ckl = c.c[0][k]*d.c[1][l];
		Csp[kl] = d.c[0][l]/d.c[1][l];
	    } else {
		ckl = c.c[0][k]*d.c[0][l];
	    }
	    for(int j = 0, ij = 0; j < b.K; ++j) {
		for(int i = 0; i < a.K; ++i, ++ij, ++ijkl) {
		    double A = braket[ij].x+braket[ij].y;
		    //rho doesn't include A or B multiplication
		    AB2rho[ijkl].y = 1.0/(A + B);
		    double C = a.c[0][i]*b.c[0][j]*ckl;
		    AB2rho[ijkl].x = (1.0/(sqrt(A+B)*A*B))*C;
		}
	    }
	}
    } 

    double rbra[6], rket[6];

    for(int i = 0; i < 3; ++i) {
	rbra[i] = r1[i];
	rbra[i+3] = r2[i];
	rket[i] = r3[i];
	rket[i+3] = r4[i];
    }

    double *devPtr;
    // make space for I[out]
    int sizeI = 0;
    if(a.L == 2 && b.L == 0 && c.L == 0 && d.L == 0)  sizeI = 6;
    else if(a.L == 2 && b.L == 0 && c.L == 0 && d.L == 0)  sizeI = 10;
    else if(a.L == 1 && b.L == 0 && c.L == 1 && d.L == 0)  sizeI = 9;

    // size includes: exponents of bra, ket;
    // rho and AB2+contractions; Csp; center information;
    int size = Kbra*2 + Kket*2 + Kbra*Kket*2 + Kket*sp + 12 + sizeI;    
    cudaMalloc(&devPtr, size*sizeof(double));

    //put all memory onto the device
    double2 *dev_braket = (double2*)(devPtr);
    double2 *dev_AB2rho = dev_braket + Kbra + Kket;
    double *dev_Csp = (double*)(dev_AB2rho + Kbra*Kket);
    double *dev_rbra = dev_Csp + Kket*sp;
    double *dev_rket = dev_rbra + 6;
    double *dev_I = dev_rket + 6;

    //TO do: pack the values into a single continuous memory segment such that
    //only a single memory copy is performed
    cudaMemcpy(dev_braket, braket, (Kbra + Kket)*sizeof(double2), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_AB2rho, AB2rho, Kbra*Kket*sizeof(double2), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_Csp, Csp, sp*Kket*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_rbra, rbra, 6*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_rket, rket, 6*sizeof(double), cudaMemcpyHostToDevice);
    assertCudaSuccess;

    dim3 dimG = dim3(1,1);
    // map bra to y dimension and ket to z dimension s.t. threadblock is < 128
    dim3 dimB = dim3(1, Kbra, 1);
    int bsize = int(dimB.x*dimB.y*dimB.z);
    int Ns = 0;
    Ns += Kbra + Kket; //eij and ekl
    Ns += 3*Kbra + 3*Kket; // rA and rB
    Ns += Kbra; // A
    Ns += Kbra; // A1
    Ns += Kket; // B
    Ns += Kket; // B1
    if((a.L + b.L) != 0) {
	Ns += 3*Kbra; // rAi
    }
    if((c.L + d.L) != 0) {
	Ns += 3*Kket; // rBk
    }
    Ns = std::max(Ns , bsize); // make sure large enough for reduction
    if(a.L == 2 && b.L == 0 && c.L == 0 && d.L == 0) {
	cuRysq_eri_dsss<<<dimG, dimB, Ns*sizeof(double)>>>(flags, tol,
							   dev_AB2rho,
							   Kbra, Kket, dev_braket,
							   dev_rbra, dev_rket, scale, dev_I);
    }
    else if(a.L == 3 && b.L == 0 && c.L == 0 && d.L == 0) {
	cuRysq_eri_fsss<<<dimG, dimB, Ns*sizeof(double)>>>(flags, tol,
							   dev_AB2rho,
							   Kbra, Kket, dev_braket,
							   dev_rbra, dev_rket, scale, dev_I);
    }
    else if(a.L == 1 && b.L == 0 && c.L == 1 && d.L == 0) {
	cuRysq_eri_psps<<<dimG, dimB, Ns*sizeof(double)>>>(flags, tol,
							   dev_AB2rho,
							   Kbra, Kket, dev_braket,
							   dev_rbra, dev_rket, scale, dev_I);
    }
    assertCudaSuccess;

    cudaMemcpy(I, dev_I, sizeI*sizeof(double), cudaMemcpyDeviceToHost);
    assertCudaSuccess;
    cudaFree(devPtr);
    assertCudaSuccess;

    //std::cout << I[1]<< std::endl;
    return 0;
}


