#ifndef _RYSQ_CUDA_BRAKET_H_
#define _RYSQ_CUDA_BRAKET_H_

#include <math.h>
#include <cuda_runtime.h>

namespace rysq {
    namespace cuda {

	class Braket {
	public:

	    class Dimension {
	    public:
		ushort Kbra, Kket;
		unsigned char L[4];
		ushort m, n, mn, N;
		ushort mi, mj, nk, nl;
		ushort eri;
		Dimension() {}
		Dimension(const Shell::Impl::Pair &bra,
			  const Shell::Impl::Pair &ket) {
		    Kbra = bra.K();
		    Kket = ket.K();
		    
		    L[0] = bra[0].L();
		    L[1] = bra[1].L();
		    L[2] = ket[0].L();
		    L[3] = ket[1].L();

		    m = L[0] + L[1] + 1;
		    n = L[2] + L[3] + 1;
		    mn = m*n;

		    mi = L[0] + 1;
		    mj = L[1] + 1;
		    nk = L[2] + 1;
		    nl = L[3] + 1;

		    N = (bra.L() + ket.L())/2 + 1;
		    eri = bra.size()*ket.size();
		}
	    };

	    Dimension dim;
	    signed char type[4];
	    int Kbra, Kket;
	    int shuffle;
	    int hybrid;
	    size_t size;
	    char *gPtr;

	    Braket(void *gPtr) {
		this->gPtr = (char*)gPtr;
	    }

	    void set(const Shell::Impl::Pair &bra,
		     const Shell::Impl::Pair &ket,
		     int shuffle) {

		dim = Dimension(bra, ket);

		type[0] = bra[0].type;
		type[1] = bra[1].type;
		type[2] = ket[0].type;
		type[3] = ket[1].type;

		Kbra = bra.K();
		Kket = ket.K();
		size = (Kbra*Kket)*sizeof(double2);
		size += (Kbra + Kket)*sizeof(double2);
		size += Kbra*bra.hybrid()*sizeof(double);
		size += Kket*ket.hybrid()*sizeof(double);
		this->shuffle = shuffle;
		hybrid = bra.hybrid() + ket.hybrid();

		char *origin = new char[size];
		size_t bytes = 0;
		double2 *AB, *aij, *akl;
		double *Cps;

		aij = (double2*)(origin + bytes);
		bytes += setExps(bra, aij);
		akl = (double2*)(origin + bytes);
		bytes += setExps(ket, akl);
		AB = (double2*)(origin + bytes);
		bytes += setAB(bra, ket, AB);
		Cps = (double*)(origin + bytes);
		bytes += setCps(bra, Cps);
		Cps = (double*)(origin + bytes);
		bytes += setCps(ket, Cps);

		memcpy_to_device(this->gPtr, origin, bytes);
		delete origin;
	    }

// 	    ~Braket() {
// 		if (origin) cuda::free(origin);
// 		//if (gData) cudaFree(gData);
// 	    }

	    size_t setAB(const Shell::Impl::Pair &bra,
			 const Shell::Impl::Pair &ket, double2 *AB) {
		size_t L = bra.L() + ket.L();
		// ket
		for(size_t l = 0, ijkl = 0; l < ket[1].K(); ++l) {
		    for(size_t k = 0; k < ket[0].K(); ++k) {
			double B = ket[0].a[k] + ket[1].a[l];
			double Ckl = ket[0].C[0][k]*ket[1].C[0][l];
			// bra
			for(size_t j = 0; j < bra[1].K(); ++j) {
			    for(size_t i = 0; i < bra[0].K(); ++i, ++ijkl) {
				// compute values
				double A = bra[0].a[i] + bra[1].a[j];
				double Cij = bra[0].C[0][i]*bra[1].C[0][j];
				double C = Cij*Ckl;
				double AB2 = 1.0/(sqrt(A+B)*A*B);
				AB[ijkl].x = C*AB2;
				AB[ijkl].y = 1.0/(A + B);
				if (L == 0) AB[ijkl].y *= A*B;
				if (L == 1) AB[ijkl].y *= A;
			    }
			}
		    }
		}
		return bra.K()*ket.K()*sizeof(double2);
	    }

	    size_t setExps(const Shell::Impl::Pair &side, double2 *dest) {
		for(size_t j = 0, ij = 0; j < side[1].K(); ++j) {
		    for(size_t i = 0; i < side[0].K(); ++i, ++ij) {
			dest[ij].x = side[0].a[i];
			dest[ij].y = side[1].a[j];
		    }
		}
		return side.K()*sizeof(double2);
	    }

	    size_t setCps(const Shell::Impl::Pair &side, double *Cps) {
		double *Cps0 = Cps;
		double *Cps1 = Cps + (side[0].type < 0)*side[0].K();
		for(size_t j = 0, ij = 0; j < side[1].K(); ++j) {
		    double Cj = side[1].C[0][j];
		    for(size_t i = 0; i < side[0].K(); ++i, ++ij) {
			double Ci = side[0].C[0][i];
			if (side[0].type < 0) Cps0[ij] = (side[0].C[1][i])/Ci;
			if (side[1].type < 0) Cps1[ij] = (side[1].C[1][j])/Cj;
		    }
		}
		return (side.K()*side.hybrid())*sizeof(double);
	    }

	};

    }
}

#endif /* _RYSQ_CUDA_BRAKET_H_ */
