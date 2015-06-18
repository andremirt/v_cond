/**
   @file
   To do: the file has become very messy needs to be cleaned up a bit
   the device memory accesses should be converted to use texture
   there is no support for double texture fetches, however the thing can be accomplished
   using integer texture fetches and unpacking int4 into two double
   The memory allocation and copy should be simplified

   finally, the code needs to be slightly modified to have high functions on center C.
   instead of on center D.
*/

#include "eri/rysq_eri_int2d.h"
#include "cuda/curysq_util.h"
#include "vec.h"
#include "util.h"
#include <math.h>
#include "roots/rysq_roots0.h"
#include "roots/rysq_roots1.h"
#include <stdio.h>
#include <iostream>
#include <algorithm>

#define pow2(x) ((x)*(x))

extern __shared__ double shmem[];

/**
   @brief
   @param a
   @param b
*/
void cuRysq_packExp(Rysq_shell_t a, Rysq_shell_t b,
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
   @param
*/
__global__ void cuRysq_eri_ssss(int flags, double tole, 
				double2 *dev_AB,
				int Kbra, int Kket, double2 *dev_braket,
				double *dev_rbra, double *dev_rket,
				double scale, double *dev_I) {
    
    typedef struct {
	double rbra[6];
	double rket[6];
	double rij[3];
	double rkl[3];
	double rij2;
	double rkl2;
	double *eij, *ekl;
	double *rA, *rB;
    } sh_t;

    short rank = ctaRank();
    short bsize = ctaSize();

    __shared__ sh_t sh;

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
	// computed bra values
	sh.rA = shmem;
	sh.eij = sh.rA + 3*Kbra;

	// computed ket values
	sh.rB = sh.eij + Kbra;
	sh.ekl = sh.rB + 3*Kket;
    }

    __syncthreads();

    // compute bra values
    for(int Kij = rank; Kij < Kbra; Kij += bsize) {
	double A1 = 1.0/(dev_braket[Kij].x + dev_braket[Kij].y);
	sh.eij[Kij] = exp(-dev_braket[Kij].x*dev_braket[Kij].y*A1*sh.rij2);
	for(int i = 0; i<3; ++i) {
	    sh.rA[i+Kij*3] = A1*(dev_braket[Kij].x*sh.rbra[i] + dev_braket[Kij].y*sh.rbra[i+3]);
	}
    }

    // compute ket values
    for(int Kkl = rank; Kkl < Kket; Kkl += bsize) {
	double B1 = 1.0/(dev_braket[Kkl+Kbra].x + dev_braket[Kkl+Kbra].y);
	sh.ekl[Kkl] = exp(-dev_braket[Kkl+Kbra].x*dev_braket[Kkl+Kbra].y*B1*sh.rkl2);
	for(int i = 0; i<3; ++i) {
	    sh.rB[i+Kkl*3] = B1*(dev_braket[Kkl+Kbra].x*sh.rket[i] + dev_braket[Kkl+Kbra].y*sh.rket[i+3]);
	}
    }
    __syncthreads();

    double q0 = 0.0;

    // compute contractions
    for(int Kkl = threadIdx.z; Kkl < Kket; Kkl += blockDim.z) {
	//the bra Kij contractions are mapped to threadIdx.y
	int K = threadIdx.y + Kkl*Kbra;

	double2 AB = dev_AB[K];
	double CAB2e = AB.x*sh.eij[threadIdx.y]*sh.ekl[Kkl];

	//if (fabs(CAB2e) < tole) continue;
	
	double X = 0.0;
	for(int i = 0; i < 3; ++i) {
	    double rAB = sh.rA[i+threadIdx.y*3] - sh.rB[i+Kkl*3];
	    X += pow2(rAB);
	}		
	X *= AB.y;

	double W = Rysq_roots0(X);
	q0 += CAB2e*W;

    }
    __syncthreads();
    
    // reduce values
    q0 = ctaReduce(q0,  shmem);

    // put values into shared memory
    if(rank == 0) {
	// scale and put values into global memory
	dev_I[0] = scale*SQRT_4PI5*q0;
    }
}

/**
   @brief
   @param flags
*/
__global__ void cuRysq_eri_sssp(int flags, double tole, 
				double2 *dev_AB2rho,
				int Kbra, int Kket, double2 *dev_braket,
				double *dev_rbra, double *dev_rket,
				double scale, double *dev_I) {
    
    typedef struct {
	double rbra[6];
	double rket[6];
	double rij[3];
	double rkl[3];
	double rij2;
	double rkl2;
	double *eij, *ekl;
	double *rA, *rB;
	double *B;
	double *rBk;
    } sh_t;

    short rank = ctaRank();
    short bsize = ctaSize();

    __shared__ sh_t sh;

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
	// computed bra values
	sh.rA = shmem;
	sh.eij = sh.rA + 3*Kbra;

	// computed ket values
	sh.rB = sh.eij + Kbra;
	sh.ekl = sh.rB + 3*Kket;

	sh.B = sh.ekl + Kket;
	sh.rBk = sh.B + Kket;
    }
    __syncthreads();

    // compute bra values
    for(int Kij = rank; Kij < Kbra; Kij += bsize) {
	double A1 = 1.0/(dev_braket[Kij].x + dev_braket[Kij].y);
	sh.eij[Kij] = exp(-dev_braket[Kij].x*dev_braket[Kij].y*A1*sh.rij2);
	for(int i=0; i<3; ++i) {
	    sh.rA[i+Kij*3] = A1*(dev_braket[Kij].x*sh.rbra[i] + dev_braket[Kij].y*sh.rbra[i+3]);
	}
    }

    // compute ket values
    for(int Kkl = rank; Kkl < Kket; Kkl += bsize) {
	sh.B[Kkl] = dev_braket[Kkl+Kbra].x + dev_braket[Kkl+Kbra].y;
	double B1 = 1.0/(dev_braket[Kkl+Kbra].x + dev_braket[Kkl+Kbra].y);
	sh.ekl[Kkl] = exp(-dev_braket[Kkl+Kbra].x*dev_braket[Kkl+Kbra].y*B1*sh.rkl2);
	for(int i = 0; i < 3; ++i) {
	    sh.rB[i+Kkl*3] = B1*(dev_braket[Kkl+Kbra].x*sh.rket[i] + dev_braket[Kkl+Kbra].y*sh.rket[i+3]);
	    sh.rBk[i+Kkl*3] = sh.rB[i+Kkl*3] - sh.rket[i+3];
	}
    }
    __syncthreads();

    double q[3] = { 0.0, 0.0, 0.0 };

    // compute contractions
    for(unsigned short int Kkl = threadIdx.z; Kkl < Kket; Kkl += blockDim.z) {
	//the bra Kij contractions are mapped to threadIdx.y
	const unsigned short Kij = threadIdx.y;
	const unsigned short K = Kij + Kkl*Kbra;

	double CAB2e = dev_AB2rho[K].x*sh.eij[Kij]*sh.ekl[Kkl];
	// should be absolute value to compare
	//if (fabs(CAB2e) < tole) continue;
	
	double X = 0.0;
	for(int i = 0; i < 3; ++i) {
	    double rAB = (sh.rA[i+Kij*3]-sh.rB[i+Kkl*3]);
	    X += pow2(rAB);
	}		
	X *= (sh.B[Kkl])*dev_AB2rho[K].y;
	double t2, W;
	Rysq_roots1(X, &t2, &W);

	CAB2e *= W;
	t2 *= dev_AB2rho[K].y;
	
	for(int i = 0; i < 3; ++i) {
	    q[i] += CAB2e*(sh.rBk[i+Kkl*3] + (sh.rA[i+Kij*3]-sh.rB[i+Kkl*3])*t2);
	}
    }
    __syncthreads();

    // reduce values
    q[0] = ctaReduce(q[0],  shmem);
    q[1] = ctaReduce(q[1],  shmem);
    q[2] = ctaReduce(q[2],  shmem);

    // put values into shared memory
    if(rank == 0) {
	for(int i = 0; i < 3; ++i) {
	    shmem[i] = q[i];
	}
    }
    __syncthreads();

    // scale and put values into global memory
    for(int i = rank; i < 3; i += bsize) {
	dev_I[i] = scale*SQRT_4PI5*shmem[i];
    }
}


/**
   @brief
   @param flags
   @param tole
   @param bra
   @param ket
   @param dev_rho
   @param dev_AB2
   @param dev_rbra
   @param dev_rket
   @param scale
   @param dev_I
*/
__global__ void cuRysq_eri_ssssp(int flags, double tole, 
				 double2 *dev_AB2rho, double *dev_Csp,
				 int Kbra, int Kket, double2 *dev_braket,
				 double *dev_rbra, double *dev_rket,
				 double scale, double *dev_I) {

    
    typedef struct {
	double rbra[6];
	double rket[6];
	double rij[3];
	double rkl[3];
	double rij2;
	double rkl2;
	double *eij, *ekl;
	double *rA, *rB;
	double *B;
	double *rBk;
    } sh_t;

    short rank = ctaRank();
    short bsize = ctaSize();

    __shared__ sh_t sh;

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
	// computed bra values
	sh.rA = shmem;
	sh.eij = sh.rA + 3*Kbra;

	// computed ket values
	sh.rB = sh.eij + Kbra;
	sh.ekl = sh.rB + 3*Kket;

	sh.B = sh.ekl + Kket;
	sh.rBk = sh.B + Kket;
    }
    __syncthreads();

    // compute bra values
    for(int Kij = rank; Kij < Kbra; Kij += bsize) {
	double A1 = 1.0/(dev_braket[Kij].x + dev_braket[Kij].y);
	sh.eij[Kij] = exp(-dev_braket[Kij].x*dev_braket[Kij].y*A1*sh.rij2);
	for(int i=0; i<3; ++i) {
	    sh.rA[i+Kij*3] = A1*(dev_braket[Kij].x*sh.rbra[i] + dev_braket[Kij].y*sh.rbra[i+3]);
	}
    }

    // compute ket values
    for(int Kkl = rank; Kkl < Kket; Kkl += bsize) {
	sh.B[Kkl] = dev_braket[Kkl+Kbra].x + dev_braket[Kkl+Kbra].y;
	double B1 = 1.0/(dev_braket[Kkl+Kbra].x + dev_braket[Kkl+Kbra].y);
	sh.ekl[Kkl] = exp(-dev_braket[Kkl+Kbra].x*dev_braket[Kkl+Kbra].y*B1*sh.rkl2);
	for(int i = 0; i < 3; ++i) {
	    sh.rB[i+Kkl*3] = B1*(dev_braket[Kkl+Kbra].x*sh.rket[i] + dev_braket[Kkl+Kbra].y*sh.rket[i+3]);
	    sh.rBk[i+Kkl*3] = sh.rB[i+Kkl*3] - sh.rket[i+3];
	}
    }
    __syncthreads();

    double q[4] = { 0.0, 0.0, 0.0, 0.0 };

    // compute contractions
    for(unsigned short int Kkl = threadIdx.z; Kkl < Kket; Kkl += blockDim.z) {
	//the bra Kij contractions are mapped to threadIdx.y
	const unsigned short Kij = threadIdx.y;
	const unsigned short K = Kij + Kkl*Kbra;
       
	double CAB2e = dev_AB2rho[K].x*sh.eij[Kij]*sh.ekl[Kkl];

	//if (fabs(CAB2e) < tole) continue;
       
	double X = 0.0;
	for(int i = 0; i < 3; ++i) {
	    double rAB = (sh.rA[i+Kij*3]-sh.rB[i+Kkl*3]);
	    X += pow2(rAB);
	}		
	X *= (sh.B[Kkl])*dev_AB2rho[K].y;
	double t2, W;
	Rysq_roots1(X, &t2, &W);

	CAB2e *= W;
	t2 *= dev_AB2rho[K].y;

	q[0] += CAB2e*dev_Csp[Kkl];

	for(int i = 0; i < 3; ++i) {
	    q[i+1] += CAB2e*(sh.rBk[i+Kkl*3] + (sh.rA[i+Kij*3]-sh.rB[i+Kkl*3])*t2);
	}
    }
    __syncthreads();

    // reduce values
    q[0] = ctaReduce(q[0],  shmem);
    q[1] = ctaReduce(q[1],  shmem);
    q[2] = ctaReduce(q[2],  shmem);
    q[3] = ctaReduce(q[3],  shmem);

    // put values into shared memory
    if(rank == 0) {
	for(int i = 0; i < 4; ++i) {
	    shmem[i] = q[i];
	}
    }
    __syncthreads();

    // scale and put values into global memory
    for(int i = rank; i < 4; i += bsize) {
	dev_I[i] = scale*SQRT_4PI5*shmem[i];
    }
}

/**
   @brief
*/
int cuRysq_eri1(int flags, double tol,
		Rysq_shell_t a, double *ri,
		Rysq_shell_t b, double *rj,
		Rysq_shell_t c, double *rk,
		Rysq_shell_t d, double *rl,
		double scale, double *I) {

    int L = a.L + b.L + c.L + d.L;
    int mask = mask(0,1,2,3);
    if(L == 1) {
	if(a.L == 1)  mask = mask(2,3,1,0);
	else if(b.L == 1) mask = mask(2,3,0,1);
	else if(c.L == 1) mask = mask(0,1,3,2);
    }

    double *r1 = ri; 
    double *r2 = rj;
    double *r3 = rk; 
    double *r4 = rl;
    shuffle(a,b,c,d,mask);
    shuffle(r1,r2,r3,r4,mask);

    bool sp = (d.nc == 2);

    int Kbra = a.K*b.K;
    int Kket = c.K*d.K;

    double2 braket[Kbra+Kket];
    cuRysq_packExp(a, b, c, d, braket);

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
		    //rho doesn't include B multiplication
		    AB2rho[ijkl].y = A/(A + B);
		    if(L==0) AB2rho[ijkl].y *= B;
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
    if(L == 0) sizeI = 1; 
    else if(sp) sizeI = 4;
    else sizeI = 3;
    // size includes: exponents of bra, ket;
    // rho and AB2+contractions; Csp; center information;
    int size = Kbra*2 + Kket*2 + Kbra*Kket*2 + Kket*sp + 12 + sizeI;    
    cudaMalloc(&devPtr, size*sizeof(double));
    assertCudaSuccess;

    //put all memory onto the device
    double2 *dev_braket = (double2*)(devPtr);
    double2 *dev_AB2rho = dev_braket + Kbra + Kket;
    double *dev_Csp = (double*)(dev_AB2rho + Kbra*Kket);
    double *dev_rbra = dev_Csp + Kket*sp;
    double *dev_rket = dev_rbra + 6;
    double *dev_I = dev_rket + 6;

    cudaMemcpy(dev_braket, braket, (Kbra + Kket)*sizeof(double2), cudaMemcpyHostToDevice);
    assertCudaSuccess;
    cudaMemcpy(dev_AB2rho, AB2rho, Kbra*Kket*sizeof(double2), cudaMemcpyHostToDevice);
    assertCudaSuccess;
    cudaMemcpy(dev_Csp, Csp, sp*Kket*sizeof(double), cudaMemcpyHostToDevice);
    assertCudaSuccess;
    cudaMemcpy(dev_rbra, rbra, 6*sizeof(double), cudaMemcpyHostToDevice);
    assertCudaSuccess;
    cudaMemcpy(dev_rket, rket, 6*sizeof(double), cudaMemcpyHostToDevice);
    assertCudaSuccess;

    dim3 dimG = dim3(1,1);
    // map bra to y dimension and ket to z dimension s.t. threadblock is < 128
    dim3 dimB = dim3(1, Kbra, 1);
    int Ns = 0;
    Ns += Kbra + Kket; //eij and ekl
    Ns += 3*Kbra + 3*Kket; // rA and rB
    if(L == 1) {
	Ns += Kket; // B
	Ns += 3*Kket; // rBk
    }
    Ns = std::max(Ns , int(dimB.x*dimB.y*dimB.z)); // make sure large enough for reduction
    if(L == 0) { // <ss|ss>
	cuRysq_eri_ssss<<<dimG, dimB, Ns*sizeof(double)>>>(flags, tol,
							   dev_AB2rho, Kbra, Kket, dev_braket,
							   dev_rbra, dev_rket, scale, dev_I);
    }
    else if(sp) { // <ss|ssp>
	cuRysq_eri_ssssp<<<dimG, dimB, Ns*sizeof(double)>>>(flags, tol,
							    dev_AB2rho, dev_Csp,
							    Kbra, Kket, dev_braket,
							    dev_rbra, dev_rket, scale, dev_I);
    } else { // <ss|sp>
	cuRysq_eri_sssp<<<dimG, dimB, Ns*sizeof(double)>>>(flags, tol,
							   dev_AB2rho, Kbra, Kket, dev_braket,
							   dev_rbra, dev_rket, scale, dev_I);
    }
    assertCudaSuccess;
    
    cudaMemcpy(I, dev_I, sizeI*sizeof(double), cudaMemcpyDeviceToHost);
    assertCudaSuccess;
    cudaFree(devPtr);
    assertCudaSuccess;

    return 0;
}

