// #include "rysq.hpp"
// #include "cuda/host.hpp"
// #include "cuda/kernels.hpp"
// #include "cuda/kernels/side.hpp"
// #include "cuda/kernels/device.hpp"
// #include "roots/rysq_roots.h"
// #include <stdio.h>

// using namespace rysq;
// using namespace rysq::cuda::kernels;

// __constant__  double2 _psss_cBraket[36*36 + 2*36 + 36];

// extern __shared__ double dmem[];

// __global__ 
// static void psss(ushort Kbra, ushort Kket,
// 		 const rysq::array::adapter<double,3> gCenters,
// 		 const rysq::array::adapter<int,4> gQuartets, ushort shuffle,
// 		 double *gEri, double cutoff);

// __global__ 
// static void spsss(ushort Kbra, ushort Kket,
// 		 const rysq::array::adapter<double,3> gCenters,
// 		 const rysq::array::adapter<int,4> gQuartets, ushort shuffle,
// 		 double *gEri, double cutoff);

// void rysq::cuda::kernels::psss(const Braket &data,
// 			       const rysq::array::adapter<double,3> gCenters,
// 			       const rysq::array::adapter<int,4> gQuartets,
// 			       double *gEri, const Parameters &parameters) {

//     ushort Kbra = data.Kbra;
//     ushort Kket = data.Kket;
// //     ushort K = Kbra*Kket;
//     dim3 gDim = dim3(gQuartets.length);
//     dim3 bDim = dim3(1, Kbra, Kket); // map bra to y dimension
//     bDim.z = std::min(bDim.z, 128/bDim.y); // maximum block size
//     size_t bSize = bDim.x*bDim.y*bDim.z;

//     cudaMemcpyToSymbol("_psss_cBraket", data.gPtr, data.size, 0,
// 		       cudaMemcpyDeviceToDevice);
//     assert_cuda_success(0);

//     size_t Ns = 0;
//     Ns += Kbra + Kket; //eij and ekl
//     Ns += 3*(Kbra + Kket); // rA and rB
//     Ns += 4*Kket; // rBk + B
//     Ns += data.hybrid*Kket; // Cp/Cs coefficient

//     Ns = std::max(Ns, bSize/2 + bSize%2); // make sure large enough for reduction
//     Ns *= sizeof(double);

//     if (data.hybrid) {
// 	::spsss<<<gDim, bDim, Ns>>>(Kbra, Kket, gCenters, gQuartets,
// 				    data.shuffle, gEri, 0);
//     }
//     else {
// 	::psss<<<gDim, bDim, Ns>>>(Kbra, Kket, gCenters, gQuartets,
// 				   data.shuffle, gEri, 0);
//     }

//     cudaThreadSynchronize();
//     assert_cuda_success(0);

// }


// __global__ 
// static void psss(ushort Kbra, ushort Kket,
// 		 const rysq::array::adapter<double,3> gCenters,
// 		 const rysq::array::adapter<int,4> gQuartets, ushort shuffle,
// 		 double *gEri, double cutoff) {

//     using namespace device;

//     const int SIZE = 3;

//     ushort bRank = block::rank();
//     ushort bSize = block::size();

//     __shared__ Side<0,0> bra;
//     __shared__ Side<1,0> ket;

//     if(bRank == 0)  {
// 	size_t offset = initialize(Kbra, dmem, bra);
// 	initialize(Kket, dmem + offset, ket);
//     }

//     ushort gRank = device::grid::rank();
//     initialize(gCenters, (int4*)gQuartets[gRank], shuffle,
// 	       bra, ket, bRank, bSize);
//     __syncthreads();

//     {
// 	double2 *cAij = _psss_cBraket;
// 	double2 *cAkl = cAij + bra.K;
// 	initialize(cAij, bra, bRank, bSize);
// 	initialize(cAkl, ket, bRank, bSize);
//     }
//     double2 *cAB = _psss_cBraket + (bra.K + ket.K);

//     __syncthreads();

//     double q[3] = { 0.0, 0.0, 0.0 };

//     // compute contractions
//     for (ushort Kkl = threadIdx.z; Kkl < ket.K; Kkl += blockDim.z) {
// 	//the bra Kij contractions are mapped to threadIdx.y
// 	int Kij = threadIdx.y;
// 	int K =  Kij + Kkl*bra.K;

// 	double2 AB = cAB[K];
// 	double CAB2e = AB.x*bra.e[Kij]*ket.e[Kkl];

// 	if (fabs(CAB2e) < cutoff) continue;
	
// 	double X = distance2(&bra.rA[Kij*3], &ket.rA[Kkl*3]);
// 	X *= (ket.A[Kkl])*AB.y;

// 	double t2, W;
// 	rysq::roots<1>(X, &t2, &W);

// 	CAB2e *= W;
// 	t2 *= AB.y;
// 	for(int i = 0; i < 3; ++i) {
// 	    q[i] += CAB2e*(ket.rAi[i+Kkl*3] + t2*(bra.rA[i+Kij*3] - ket.rA[i+Kkl*3]));
// 	}

//     }
//     __syncthreads();
    
//     // reduce values
//     q[0] = block::reduce(q[0], dmem);
//     q[1] = block::reduce(q[1], dmem);
//     q[2] = block::reduce(q[2], dmem);

//     __shared__ double Q[SIZE];
//     if (bRank == 0) copy(Q, q, SIZE);
//     scale(gEri + grid::rank()*SIZE, SQRT_4PI5, Q, SIZE, bSize, bRank);

// }

// __global__ 
// static void spsss(ushort Kbra, ushort Kket,
// 		 const rysq::array::adapter<double,3> gCenters,
// 		 const rysq::array::adapter<int,4> gQuartets, ushort shuffle,
// 		 double *gEri, double cutoff) {

//     using namespace device;

//     const int SIZE = 4;

//     ushort bRank = block::rank();
//     ushort bSize = block::size();
    
//     __shared__ Side<0,0> bra;
//     __shared__ Side<1,0> ket;
//     __shared__ double *Cps;

//     if(bRank == 0)  {
// 	size_t offset = initialize(Kbra, dmem, bra);
// 	offset += initialize(Kket, dmem + offset, ket);
// 	Cps = dmem + offset;
//     }

//     ushort gRank = device::grid::rank();
//     initialize(gCenters, (int4*)gQuartets[gRank], shuffle,
// 	       bra, ket, bRank, bSize);
//     __syncthreads();

//     {
// 	double2 *cAij = _psss_cBraket;
// 	initialize(cAij, bra, bRank, bSize);
// 	double2 *cAkl = _psss_cBraket + bra.K;
// 	initialize(cAkl, ket, bRank, bSize);
//     }
//     double2 *cAB = _psss_cBraket + (bra.K + ket.K);

//     double *cCps = (double*)(cAB + bra.K*ket.K);
//     copy(Cps, cCps, ket.K, bSize, bRank);

//     __syncthreads();

//     double q[SIZE] = { 0.0, 0.0, 0.0, 0.0 };

//     // compute contractions
//     for (ushort Kkl = threadIdx.z; Kkl < ket.K; Kkl += blockDim.z) {
// 	//the bra Kij contractions are mapped to threadIdx.y
// 	int Kij = threadIdx.y;
// 	int K =  Kij + Kkl*bra.K;

// 	double2 AB = cAB[K];
// 	double CAB2e = AB.x*bra.e[Kij]*ket.e[Kkl];
// 	double C = Cps[Kkl];
// 	if (max(fabs(CAB2e), fabs(C*CAB2e)) < cutoff) continue;
	
// 	double X = distance2(&bra.rA[Kij*3], &ket.rA[Kkl*3]);
// 	X *= (ket.A[Kkl])*AB.y;

// 	double t2, W;
// 	rysq::roots<1>(X, &t2, &W);

// 	CAB2e *= W;
// 	t2 *= AB.y;

// 	q[0] += CAB2e;
// 	CAB2e *= C;
// 	for(int i = 0; i < 3; ++i) {
// 	    q[i+1] += CAB2e*(ket.rAi[i+Kkl*3] + t2*(bra.rA[i+Kij*3] - ket.rA[i+Kkl*3]));
// 	}

//     }
//     __syncthreads();
    
//     // reduce values
//     q[0] = block::reduce(q[0], dmem);
//     q[1] = block::reduce(q[1], dmem);
//     q[2] = block::reduce(q[2], dmem);
//     q[3] = block::reduce(q[3], dmem);

//     __shared__ double Q[SIZE];
//     if (bRank == 0) copy(Q, q, SIZE);
//     scale(gEri + grid::rank()*SIZE, SQRT_4PI5, Q, SIZE, bSize, bRank);

// }
