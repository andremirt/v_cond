// #include "rysq.hpp"
// #include "cuda/host.hpp"
// #include "cuda/kernels.hpp"
// #include "cuda/kernels/side.hpp"
// #include "cuda/kernels/device.hpp"
// #include "roots/rysq_roots.h"

// using namespace rysq::cuda::kernels;

// __device__ __constant__ double2 _ssss_cBraket[36*36 + 2*36];

// extern __shared__ double dmem[];

// __global__ 
// static void ssss(ushort Kbra, ushort Kket,
// 		 const rysq::array::adapter<double,3> gCenters,
// 		 const rysq::array::adapter<int,4> gQuartets,
// 		 ushort shuffle, double *gEri, double cutoff);

// using namespace  rysq;

// void rysq::cuda::kernels::ssss(const Braket &data,
// 			       const rysq::array::adapter<double,3> gCenters,
// 			       const rysq::array::adapter<int,4> gQuartets,
// 			       double *gEri, const Parameters &parameters) {

//     ushort Kbra = data.Kbra;
//     ushort Kket = data.Kket;
//     //ushort K = Kbra*Kket;
//     dim3 gDim = dim3(gQuartets.length);
//     dim3 bDim = dim3(1, Kbra, Kket); // map bra to y dimension
//     bDim.z = std::min(bDim.z, 128/bDim.y); // maximum block size
//     size_t bSize = bDim.x*bDim.y*bDim.z;

//     cudaMemcpyToSymbol("_ssss_cBraket", data.gPtr, data.size, 0,
// 		       cudaMemcpyDeviceToDevice);
//     assert_cuda_success(0);

//     size_t Ns = 0;
//     Ns += Kbra + Kket; //eij and ekl
//     Ns += 3*(Kbra + Kket); // rA and rB
//     Ns = std::max(Ns, bSize/2 + bSize%2); // make sure large enough for reduction
//     Ns *= sizeof(double);

//     ::ssss<<<gDim, bDim, Ns>>>(Kbra, Kket, gCenters, gQuartets,
// 			       data.shuffle, gEri, 0);
//     cudaThreadSynchronize();
//     assert_cuda_success(0);

// }


// __global__ 
// static void ssss(ushort Kbra, ushort Kket,
// 		 const rysq::array::adapter<double,3> gCenters,
// 		 const rysq::array::adapter<int,4> gQuartets,
// 		 ushort shuffle,
// 		 double *gEri, double cutoff) {

//     using namespace device;

//     ushort bRank = block::rank();
//     ushort bSize = block::size();

//     __shared__ Side<0,0> bra, ket;

//     if(bRank == 0)  {
// 	size_t offset = initialize(Kbra, dmem, bra);
// 	initialize(Kket, dmem + offset, ket);
//     }

//     ushort gRank = device::grid::rank();
//     initialize(gCenters, (int4*)gQuartets[gRank], shuffle, bra, ket, bRank, bSize);
//     __syncthreads();

//     {
// 	double2 *cAij = _ssss_cBraket;
// 	double2 *cAkl = cAij + bra.K;
// 	initialize(cAij, bra, bRank, bSize);
// 	initialize(cAkl, ket, bRank, bSize);
//     }
//     double2 *cAB = _ssss_cBraket + (bra.K + ket.K);

//     __syncthreads();

//     double q0 = 0.0;

//     // compute contractions
//     for (ushort Kkl = threadIdx.z; Kkl < ket.K; Kkl += blockDim.z) {
// 	//the bra Kij contractions are mapped to threadIdx.y
// 	int Kij = threadIdx.y;
// 	int K =  Kij + Kkl*bra.K;

// 	double2 AB = cAB[K];
// 	double CAB2e = AB.x*bra.e[Kij]*ket.e[Kkl];

// 	if (fabs(CAB2e) < cutoff) continue;
	
// 	double X = distance2(&bra.rA[Kij*3], &ket.rA[Kkl*3]);
// 	X *= AB.y;

// 	double W;
// 	rysq::roots<0>(X, NULL, &W);
// 	q0 += CAB2e*W;

//     }
//     __syncthreads();
    
//     // reduce values
//     q0 = block::reduce(q0, dmem);

//     // put values into smemared memory
//     if(bRank == 0) gEri[grid::rank()] = SQRT_4PI5*q0;
    
// }


