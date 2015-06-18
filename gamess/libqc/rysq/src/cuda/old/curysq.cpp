/**
   @file
   @brief Public function implemention for CUDA Rysq
 */

#include "rysq.h"
#include "cuda/curysq_util.h"
#include "curysq_eri.h"
#include "curysq_jk.h"
#include "vec.h"
#include "util.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>


/** @see rysq.h */
int cuRysq_init(int device) {
    cudaSetDevice(device);
    assertCudaSuccess;
    return 0;
}

/** @see rysq.h */
int cuRysq_eri(int flags, double tol,
	       Rysq_shell_t a, double *ri,
	       Rysq_shell_t b, double *rj,
	       Rysq_shell_t c, double *rk,
	       Rysq_shell_t d, double *rl,
	       int n, int index4[], double scale, double *I) {

    int dev;
    struct cudaDeviceProp prop;
    assert(cudaGetDevice(&dev) == cudaSuccess);
    assert(cudaGetDeviceProperties(&prop,dev) == cudaSuccess);

    int size2D = (a.L+1)*(b.L+1)*(c.L+1)*(d.L+1);
    int N = Rysq_num_roots(a.L, b.L, c.L, d.L);

    // allocate device memory
    double *dev_mem;
    size_t mem = 0;
    mem += N*size2D*n*3;
    mem += Rysq_eri_size(a, b, c, d)*n;
    mem *= sizeof(double);
    cudaMalloc((void**)&dev_mem, mem);
    assertCudaSuccess;

    // offsets 
    double *dev_Ix = dev_mem;
    double *dev_Iy = dev_Ix + N*size2D*n;
    double *dev_Iz = dev_Iy + N*size2D*n;
    double *dev_I = dev_Iz + N*size2D*n;

//     double *rij = (double*)malloc((nij+nkl)*6*sizeof(double));
//     double *rkl = &rij[nij*6];

//     for(int j = 0, ij = 0; j < nj; ++j) {
// 	for(int i = 0; i < ni; ++i, ij += 6) {
// 	    Vec_copy3(&ri[i*3], &rij[ij+0]);
// 	    Vec_copy3(&rj[j*3], &rij[ij+3]);
// 	}
//     }

//     for(int l = 0, kl = 0; l < nl; ++l) {
// 	for(int k = 0; k < nk; ++k, kl += 6) {
// 	    Vec_copy3(&rk[k*3], &rkl[kl+0]);
// 	    Vec_copy3(&rl[l*3], &rkl[kl+3]);
// 	}
//     }

//     cudaMemcpy(dev_rij, rij, (nij+nkl)*6*sizeof(double), cudaMemcpyHostToDevice);
//     assertCudaSuccess;    
    
    int block2 = Rysq_shell_size(a)*Rysq_shell_size(b);
    int block3 = block2*Rysq_shell_size(c);
    int block4 = block3*Rysq_shell_size(d);
    
//     if(block4 <= prop.maxThreadsPerBlock) {
// 	cuRysq_eri4(flags, tol, a, b, c, d, n,
// 		    dev_Ix, dev_Iy, dev_Iz, scale, dev_I);
//     }
//     else if(block3 <= prop.maxThreadsPerBlock) {
// 	cuRysq_eri3(flags, tol, a, b, c, d, n,
// 		    dev_Ix, dev_Iy, dev_Iz, scale, dev_I);
//     }
//     else {
// 	cuRysq_eri2(flags, tol, a, b, c, d, n,
// 		    dev_Ix, dev_Iy, dev_Iz, scale, dev_I);
//     }

    // sync host and device
    cudaThreadSynchronize();
    assertCudaSuccess;
    
    cudaMemcpy(I, dev_I, Rysq_eri_size(a, b, c, d)*n*sizeof(double), 
	       cudaMemcpyDeviceToHost);
    assertCudaSuccess;

    // free memory
    //    free(rij);
    cudaFree(dev_mem);
    assertCudaSuccess; 

    return 0;
}

/** @see rysq.h */
int cuRysq_jk(int flags, double tol,
	      Rysq_shell_t a, int ni, double *ri,
	      Rysq_shell_t b, int nj, double *rj,
	      Rysq_shell_t c, int nk, double *rk,
	      Rysq_shell_t d, int nl, double *rl,
	      double *Dij, double *Dkl, double *Dil, double *Djk,
	      double alpha,
	      double *Fij, double *Fkl, double *Fil, double *Fjk) {

    int dev;
    struct cudaDeviceProp prop;
    assert(cudaGetDevice(&dev) == cudaSuccess);
    assert(cudaGetDeviceProperties(&prop,dev) == cudaSuccess);

    int size2D = (a.L+1)*(b.L+1)*(c.L+1)*(d.L+1);
    int nij = ni*nj;
    int nkl = nk*nl;
    int N = Rysq_num_roots(a.L, b.L, c.L, d.L);

    // allocate device memory
    double *dev_mem;
    size_t mem = 0;
    mem += (nij + nkl)*6; //rij and rkl
    mem += N*size2D*nij*nkl*3; //Ix, Iy, and Iz
    mem += 2*a.size*b.size*nij; //Dij and Fij
    mem += 2*c.size*d.size*nkl; //Dkl and Fkl
    mem += 2*a.size*d.size*ni*nl; //Dil and Fil
    mem += 2*b.size*c.size*nj*nk; //Djk and Fjk
    mem *= sizeof(double);
    cudaMalloc((void**)&dev_mem, mem);
    assertCudaSuccess;

    // offsets 
    double *dev_rij = dev_mem;
    double *dev_rkl = dev_rij + nij*6;
    double *dev_Ix = dev_rkl + nkl*6;
    double *dev_Iy = dev_Ix + N*size2D*nij*nkl;
    double *dev_Iz = dev_Iy + N*size2D*nij*nkl;
    double *dev_Dij = dev_Iz + N*size2D*nij*nkl;
    double *dev_Dkl = dev_Dij + a.size*b.size*nij;
    double *dev_Dil = dev_Dkl + c.size*d.size*nkl;
    double *dev_Djk = dev_Dil + a.size*d.size*ni*nl;
    double *dev_Fij = dev_Djk + b.size*c.size*nj*nk;
    double *dev_Fkl = dev_Fij + a.size*b.size*nij;
    double *dev_Fil = dev_Fkl + c.size*d.size*nkl;
    double *dev_Fjk = dev_Fil + a.size*d.size*ni*nl;

    double *rij = (double*)malloc((nij+nkl)*6*sizeof(double));
    double *rkl = &rij[nij*6];

    for(int j = 0, ij = 0; j < nj; ++j) {
	for(int i = 0; i < ni; ++i, ij += 6) {
	    Vec_copy3(&ri[i*3], &rij[ij+0]);
	    Vec_copy3(&rj[j*3], &rij[ij+3]);
	}
    }

    for(int l = 0, kl = 0; l < nl; ++l) {
	for(int k = 0; k < nk; ++k, kl += 6) {
	    Vec_copy3(&rk[k*3], &rkl[kl+0]);
	    Vec_copy3(&rl[l*3], &rkl[kl+3]);
	}
    }

    cudaMemcpy(dev_rij, rij, (nij+nkl)*6*sizeof(double), cudaMemcpyHostToDevice);
    assertCudaSuccess;
    
    cudaMemcpy(dev_Dij, Dij, a.size*b.size*nij*sizeof(double), cudaMemcpyHostToDevice);
    assertCudaSuccess;

    cudaMemcpy(dev_Dkl, Dkl, c.size*d.size*nkl*sizeof(double), cudaMemcpyHostToDevice);
    assertCudaSuccess;

    cudaMemcpy(dev_Dil, Dil, a.size*d.size*ni*nl*sizeof(double), cudaMemcpyHostToDevice);
    assertCudaSuccess;

    cudaMemcpy(dev_Djk, Djk, b.size*c.size*nj*nk*sizeof(double), cudaMemcpyHostToDevice);
    assertCudaSuccess;

//     int block2 = Rysq_shell_size(a)*Rysq_shell_size(b);
//     int block3 = block2*Rysq_shell_size(c);
//     int block4 = block3*Rysq_shell_size(d);
    
//     if(block4 <= prop.maxThreadsPerBlock) {
// 	cuRysq_eri4(flags, tol, a, b, ni, nj, dev_rij, c, d, nk, nl, dev_rkl,
// 		    dev_Ix, dev_Iy, dev_Iz, alpha, dev_I);
//     }
//     else if(block3 <= prop.maxThreadsPerBlock) {
// 	cuRysq_eri3(flags, tol, a, b, ni, nj, dev_rij, c, d, nk, nl, dev_rkl,
// 		    dev_Ix, dev_Iy, dev_Iz, alpha, dev_I);
//     }
//     else {
// 	cuRysq_eri2(flags, tol, a, b, ni, nj, dev_rij, c, d, nk, nl, dev_rkl,
// 		    dev_Ix, dev_Iy, dev_Iz, alpha, dev_I);
//     }

//     cuRysq_jk4(flags, tol, a, b, ni, nj, dev_rij, c, d, nk, nl, dev_rkl,
// 	       dev_Ix, dev_Iy, dev_Iz, dev_Dij, dev_Dkl, dev_Dil, dev_Djk,
// 	       alpha, dev_Fij, dev_Fkl, dev_Fil, dev_Fjk);
//     assertCudaSuccess;

    // sync host and device
    cudaThreadSynchronize();
    assertCudaSuccess;
    
    cudaMemcpy(Fij, dev_Fij, a.size*b.size*nij*sizeof(double), cudaMemcpyDeviceToHost);
    assertCudaSuccess;

    cudaMemcpy(Fkl, dev_Fkl, c.size*d.size*nkl*sizeof(double), cudaMemcpyDeviceToHost);
    assertCudaSuccess;

    cudaMemcpy(Fil, dev_Fil, a.size*d.size*ni*nl*sizeof(double), cudaMemcpyDeviceToHost);
    assertCudaSuccess;

    cudaMemcpy(Fjk, dev_Fjk, b.size*c.size*nj*nk*sizeof(double), cudaMemcpyDeviceToHost);
    assertCudaSuccess;

    // free memory
    free(rij);
    cudaFree(dev_mem);
    assertCudaSuccess; 

    return 0;
}

#undef assertCudaSuccess

