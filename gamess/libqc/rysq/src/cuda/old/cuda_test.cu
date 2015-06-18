#define in(i, dim) dev_G[i + dim * inputOffset]
#define out(i, dim) dev_out[i + dim * outputOffset];

__global__ void computeG(int ni, int nj, int nk, int nl,
			 int3 constants, double *dev_G, double *dev_out) {
    short rank = ctaRank();
    short bsize = ctaSize();

    int inputOffset = std::max(ni + nj - 1, nk + nl - 1);
    int outputOffset = std::max(ni, nk) * std::max(nj, nl);

    //Loop over the 3 dimensions: x, y, z
    for (int dim = 0; dim < 3; ++dim) {

	for (int i = 0; i < std::max(ni, nk); ++i) {
	    out(i, dim) = in(i, dim);
	}
	__syncthreads();

	for (int j = 1; j < nj; ++j) {
	    for (int i = 0; i < ni + nj - 1 - j; ++i) {
		double temp = dev_G[i + 1];
		__syncthreads();
		// Figure out actual constant 'constants'
		in(i, dim) = temp + constants.dim * in(i, dim);
		__syncthreads();
		if (rank < ni) {
		    out(i + j * ni, dim) = in(i, dim);
		}
	    }
	}
    
	for (int l = 1; l < nl; ++l) {
	    for (int k = 0; k < nk + nl - 1 - l; ++k) {
		double temp = dev_G[i + 1];
		__syncthreads();
		//Figure out actual constant 'constants'
		in(k, dim) = temp + constants.dim * in(k, dim);
		__syncthreads();
		if (rank < nk) {
		    out(k + l * nk, dim) = in(k, dim);
		}
	    }
	}
    }
}
