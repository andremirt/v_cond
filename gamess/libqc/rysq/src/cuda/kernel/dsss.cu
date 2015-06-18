// #include "rysq.hpp"
// #include "cuda/host.hpp"
// #include "cuda/kernels.hpp"
// #include "cuda/kernels/side.hpp"
// #include "cuda/kernels/device.hpp"
// #include "roots/rysq_roots.h"
// #include <stdio.h>

// using namespace rysq;
// using namespace rysq::cuda::kernels;

// __constant__  double2 _dsss_cBraket[36*36 + 2*36];

// extern __shared__ double dmem[];

// __global__ 
// static void dsss(ushort Kbra, ushort Kket,
// 		 const rysq::array::adapter<double,3> gCenters,
// 		 const rysq::array::adapter<int,4> gQuartets, ushort shuffle,
// 		 double *gEri, double cutoff);

// void rysq::cuda::kernels::dsss(const Braket &data,
// 			       const rysq::array::adapter<double,3> gCenters,
// 			       const rysq::array::adapter<int,4> gQuartets,
// 			       double *gEri, const Parameters &parameters) {

//     ushort Kbra = data.Kbra;
//     ushort Kket = data.Kket;

//     dim3 gDim = dim3(gQuartets.length);
//     dim3 bDim = dim3(1, Kbra, Kket); // map bra to y dimension
//     bDim.z = std::min(bDim.z, 128/bDim.y); // maximum block size
//     size_t bSize = bDim.x*bDim.y*bDim.z;

//     cudaMemcpyToSymbol("_dsss_cBraket", data.gPtr, data.size, 0,
// 		       cudaMemcpyDeviceToDevice);
//     assert_cuda_success(0);

//     size_t Ns = 0;
//     Ns += Kbra + Kket; //eij and ekl
//     Ns += 3*(Kbra + Kket); // rA and rB
//     Ns += 5*Kbra; // rAi + A + 1/A
//     Ns += 1*Kket; // B

//     Ns = std::max(Ns, bSize/2 + bSize%2); // make sure large enough for reduction
//     Ns *= sizeof(double);

//     ::dsss<<<gDim, bDim, Ns>>>(Kbra, Kket, gCenters, gQuartets,
// 				   data.shuffle, gEri, 0);

//     cudaThreadSynchronize();
//     assert_cuda_success(0);

// }


// __global__ 
// static void dsss(ushort Kbra, ushort Kket,
// 		 const rysq::array::adapter<double,3> gCenters,
// 		 const rysq::array::adapter<int,4> gQuartets, ushort shuffle,
// 		 double *gEri, double cutoff) {

//     using namespace device;

//     const int SIZE = 6;

//     ushort bRank = block::rank();
//     ushort bSize = block::size();

//     __shared__ Side<2,0> bra;
//     __shared__ Side<0,1> ket;

//     if(bRank == 0)  {
// 	size_t offset = initialize(Kbra, dmem, bra);
// 	initialize(Kket, dmem + offset, ket);
//     }

//     ushort gRank = device::grid::rank();
//     initialize(gCenters, (int4*)gQuartets[gRank], shuffle,
// 	       bra, ket, bRank, bSize);
//     __syncthreads();

//     {
// 	double2 *cAij = _dsss_cBraket;
// 	double2 *cAkl = cAij + bra.K;
// 	initialize(cAij, bra, bRank, bSize);
// 	initialize(cAkl, ket, bRank, bSize);
//     }
//     double2 *cAB = _dsss_cBraket + (bra.K + ket.K);

//     __syncthreads();

//     double q[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

//     // compute contractions
//     for (ushort Kkl = threadIdx.z; Kkl < ket.K; Kkl += blockDim.z) {
// 	//the bra Kij contractions are mapped to threadIdx.y
// 	int Kij = threadIdx.y;
// 	int K =  Kij + Kkl*bra.K;

// 	double2 AB = cAB[K];
// 	double CAB2e = AB.x*bra.e[Kij]*ket.e[Kkl];

// 	if (fabs(CAB2e) < cutoff) continue;
	
// 	double X = distance2(&bra.rA[Kij*3], &ket.rA[Kkl*3]);
// 	//printf("X %f %f %f %f GPU\n", X, bra.A[Kij], (ket.A[Kkl]), AB.y);
// 	X *= (bra.A[Kij])*(ket.A[Kkl])*AB.y;

// 	double t2[2], W[2];
// 	rysq::roots<2>(X, t2, W);

// 	t2[0] *= AB.y;
// 	t2[1] *= AB.y;

// 	double WBt[3];
// 	WBt[0] = W[0] + W[1];
// 	WBt[1] = ket.A[Kkl]*(W[0]*t2[0] + W[1]*t2[1]);
// 	WBt[2] = pow2(ket.A[Kkl])*(W[0]*t2[0]*t2[0] + W[1]*t2[1]*t2[1]);
	
// 	double rAB = (bra.rA[Kij*3] - ket.rA[Kkl*3]);
// 	double rAi = bra.rAi[Kij*3];
// 	double t1 = 0.5*bra.A1[Kij]*(WBt[0] - WBt[1]);	

// 	//last integrals have to be swapped
// 	for(int i = 0; i < 3; ++i) {
// 	    double x1 = WBt[2]*rAB - WBt[1]*rAi;
// 	    double x2 = WBt[0]*rAi - WBt[1]*rAB;

// 	    q[i] += CAB2e*(rAB*x1 + rAi*x2 + t1);

// 	    int j = (i+1)%3;
// 	    rAB = (bra.rA[j+Kij*3] - ket.rA[j+Kkl*3]);
// 	    rAi = bra.rAi[j+Kij*3];

// 	    q[i+3] += CAB2e*(rAB*x1 + rAi*x2);
//  	}

//     }
//     __syncthreads();
    
//     // reduce values
//     q[0] = block::reduce(q[0], dmem);
//     q[1] = block::reduce(q[1], dmem);
//     q[2] = block::reduce(q[2], dmem);
//     q[3] = block::reduce(q[3], dmem);
//     q[4] = block::reduce(q[4], dmem);
//     q[5] = block::reduce(q[5], dmem);

//     __shared__ double Q[SIZE];
//     if (bRank == 0) {
// 	swap(q[4], q[5]);
// 	copy(Q, q, SIZE);
//     }
//     scale(gEri + grid::rank()*SIZE, SQRT_4PI5, Q, SIZE, bSize, bRank);

// }
