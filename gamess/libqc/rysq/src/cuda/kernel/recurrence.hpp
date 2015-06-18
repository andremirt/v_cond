#ifndef _RYSQ_CUDA_KERNEL_RECURRENCE_HPP_
#define _RYSQ_CUDA_KERNEL_RECURRENCE_HPP_

template<int N> __device__
void recurrence(const ushort &m, const ushort &n,
		const double &A, const double &B,
		const double &A1, const double &B1,
		const double *rAB, const double *rAi, const double *rBk,
	        const double *t2, const double *W, double *G,
		const ushort &bSize, const ushort &bRank) {

//     using namespace device;
//     ushort bRank = block::rank();
//     ushort bSize = block::size();

    if (bRank == 0) {
	//G[-1] = 0; // G[-1] must not be nan or inf
	ushort mn = m*n;
	for (ushort i = 0, j = 0; i < N; ++i) {
	    G[j] = 1.0;
	    j += mn;
	    G[j] = 1.0;
	    j += mn;
	    G[j] = W[i];
	    j += mn;
	}
    }

    __syncthreads();

    __shared__ double Cm[N*3], Cn[N*3];
    __shared__ double B00[N], B10[N], B01[N];

    //for (ushort a = bRank; a < N; a += bSize) {
    if (bRank < N) {
	const ushort &a = bRank;
	double t = t2[a];
	B00[a] = 0.5*t;
	B10[a] = 0.5*A1*(1.0 - B*t);
	B01[a] = 0.5*B1*(1.0 - A*t);
	for (ushort r = 0; r < 3; ++r) {
	    double q = rAB[r]*t;
	    Cm[r+a*3] = rAi[r] - B*q;
	    Cn[r+a*3] = rBk[r] + A*q;
	}
    }

    __syncthreads();

    short x = threadIdx.x;
    const ushort &y = threadIdx.y;
    const ushort &z = threadIdx.z;
    ushort yz = threadIdx.y + threadIdx.z*3;
    double *myG = G + yz*m*n;

    if (x == 0 and z < N) {
	double B = B10[z];
	double C = Cm[yz];
	myG[1] = C*myG[0];

	for(ushort i = 2; i < m; ++i) {
	    myG[i] = B*myG[i-2] + C*myG[i-1];
	    B += B10[z];
	}
    }

    __syncthreads();

    if(n == 1) return;

    if (x < m and z < N) {
	double C = Cn[yz];
	double B0 = x*B00[z];
	myG[x+m] = B0*myG[x-1] + C*myG[x];
	double q = B01[z];
	for(ushort k = 2; k < n; ++k) {
	    myG += m; // next column
	    myG[x+m] = q*myG[x-m] + B0*myG[x-1] + C*myG[x];
	    q += B01[z];
	}
    }

}


#endif /* _RYSQ_CUDA_KERNEL_RECURRENCE_HPP_ */
