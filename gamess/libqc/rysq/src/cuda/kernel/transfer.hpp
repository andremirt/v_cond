#ifndef _RYSQ_CUDA_KERNELS_TRANSFER_H_
#define _RYSQ_CUDA_KERNELS_TRANSFER_H_

// #include "cuda/kernel/device.hpp"

#include <cuda.h>
#include <cuda_runtime.h>

template<ushort N> __device__
void transfer1(const ushort &mi, const ushort &mj, const ushort &n,
	       const double *dr, double *G, double *I) {

    // using namespace device;

    const ushort &x = threadIdx.x;
    const ushort &y = threadIdx.y;
    const ushort &z = threadIdx.z;

    if (z >= N) return;

    const ushort mij = mi*mj;
    const ushort m = mi + mj - 1;

    const double dq = dr[y];
    const ushort yz = (y + z*3)*n;

    //for(ushort kl = x; kl < n; kl += blockDim.x) {
    if (x < n) {
	const ushort &kl = x;
	double *myG = G + (kl + yz)*(m);
	double *myI = I + (kl + yz)*(mij);
	// printf("%p %u %u %u\n", myI, z, y, n);
 	for(ushort i = 0; i < mi; ++i) {
	    myI[i] = myG[i];
	}
	for(ushort j = 1; j < mj; ++j) {
	    double q0 = myG[0];
	    myI += mi;
	    for(ushort i = 0; i < mi; ++i) {
		//double q0 = myG[i];
		double q1 = myG[i+1];
		q0 = dq*q0 + q1;
		myI[i] = q0;
		myG[i] = q0;
		q0 = q1;
	    }
	    for(ushort i = mi; i < (m - j); ++i) {
		//myG[i] = dq*myG[i] + myG[i+1];
		double q1 = myG[i+1];
		q0 = dq*q0 + q1;
		myG[i] = q0;
		q0 = q1;
	    }		    
	}
    }
}

template<ushort N> __device__
void transfer2(const ushort &mij, const ushort &nk, const ushort &nl,
	       const double *dr, double *G, double *I) {

    // using namespace device;

    const ushort &x = threadIdx.x;
    const ushort &y = threadIdx.y;
    const ushort &z = threadIdx.z;

    if (z >= N) return;

    const ushort nkl = nk*nl;
    const ushort n = nk + nl - 1;

    const ushort yz = (y + z*3)*mij;
    const double dq = dr[y];

    for(ushort ij = x; ij < mij; ij += blockDim.x) {
	double *myG = G + ij + yz*n;
	double *(myI) = I + ij + yz*nkl;

	for(ushort k = 0; k < nk; ++k) {
	    *(myI) = *(myG);
	    myI += mij;
	    myG += mij;
	}
	myG -= nk*mij;
	for(ushort l = 1; l < nl; ++l) {
	    double q0 = *(myG);
	    for(ushort k = 0; k < nk; ++k) {
		double q1 = *(myG + mij);
		q0 = dq*q0 + q1;
		*(myI) = q0;
		*(myG) = q0;
		q0 = q1;
		myG += mij;
		myI += mij;
	    }
	    for(ushort k = nk; k < (n - l); ++k) {
		double q1 = *(myG + mij);
		*(myG) = dq*q0 + q1;
		q0 = q1;
		myG += mij;
	    }		    
	    myG -= (n-l)*mij;
	}
    }
}


#endif /* _RYSQ_CUDA_KERNELS_TRANSFER_H_ */
