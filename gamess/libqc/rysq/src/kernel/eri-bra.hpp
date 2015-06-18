#ifndef _RYSQ_KERNEL_ERI_BRA_HPP_
#define _RYSQ_KERNEL_ERI_BRA_HPP_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <math.h>
#include <boost/mpl/int.hpp>
#include <rysq/core.hpp>
#include "vector.hpp"

#include "kernel/eri.hpp"
#include "kernel/primitives.hpp"
#include "cxx/utility/permute.hpp"

#include "kernel/transform.hpp"

BEGIN_NAMESPACE(rysq, kernel)


template<class bra_, int N>
struct eri<bra_, boost::mpl::int_<N> > : public Eri
{
    typedef bra_ bra;
    typedef void ket;
    typedef kernel::Transform<bra> Transform;
    typedef typename Transform::Data Data;
    eri(const Quartet<Shell> &quartet, Transform *transform)
	: quartet_(quartet),
	  transform_(transform),
	  primitives_(quartet) {}

    ~eri() { delete transform_; }

    void operator()(const Quartet<Center> &r,
		    Data &data,
		    const Parameters &parameters) {
	double scale = rysq::SQRT_4PI5;
	double cutoff = parameters.cutoff/(scale*quartet_.K()*10);
	apply(quartet_, r[0], r[1], r[2], r[3], scale, cutoff, primitives_,
	      (*transform_)(data));
    }

private:


    const Quartet<Shell> quartet_;
    Transform *transform_;
    eri_::Primitives<double> primitives_;

    template<typename F>
    static void apply(const Quartet<Shell> &quartet,
		      const Vector<3> &ri, const Vector<3> &rj,
		      const Vector<3> &rk, const Vector<3> &rl,
		      double scale, double cutoff,
		      eri_::Primitives<F> &primitives,
		      Transform &transform);//  {

    // 	typedef kernel::eri_::quadrature<bra> quadrature;

    // 	const Shell &a = quartet[0];
    // 	const Shell &b = quartet[1];
    // 	const Shell &c = quartet[2];
    // 	const Shell &d = quartet[3];

    // 	Vector<3> rij = ri - rj;
    // 	Vector<3> rkl = rk - rl;

    // 	double rkl2 = rkl.inner();
    // 	double rij2 = rij.inner();

    // 	double eij[a.K*b.K];
    // 	for (int Kj = 0, Kij = 0; Kj < b.K; ++Kj) {
    // 	    for (int Ki = 0; Ki < a.K; ++Ki, ++Kij) {
    // 		double A = a(Ki) + b(Kj);
    // 		double A1 = 1.0/A;
    // 		eij[Kij] = exp(-a(Ki)*b(Kj)*A1*rij2);
    // 	    }
    // 	}

    // 	int K = 0;

    // 	for (int Kl = 0; Kl < d.K; ++Kl) {
    // 	    for (int Kk = 0; Kk < c.K; ++Kk) {
    // 		double ak = c(Kk);
    // 		double al = d(Kl);
    // 		double B = ak + al;

    // 		Vector<3> rB = Vector<3>::center(ak, rk, al, rl);
    // 		Vector<3> rBk = rB - rk;

    // 		double ekl = exp(-ak*al*rkl2/B);

    // 		double Ckl[4] __attribute__ ((aligned(16)));
    // 		double Ckl_max = 0.0;
    // 		for(int l = 0, kl = 0; l < d.nc; ++l) {
    // 		    for(int k = 0; k < c.nc; ++k, ++kl) {
    // 			Ckl[kl] = c(Kk,k)*d(Kl,l);
    // 			Ckl_max = std::max(Ckl_max, fabs(Ckl[kl]));
    // 		    }
    // 		}
		    
    // 		for(int Kj = 0, Kij = 0; Kj < b.K; ++Kj) {
    // 		    for(int Ki = 0; Ki < a.K; ++Ki, ++Kij) {
    // 			double ai = a(Ki);
    // 			double aj = b(Kj); 
    // 			double A = ai + aj;
    // 			double e = eij[Ki+Kj*a.K]*ekl;
    // 			double eAB = e/(A*B*sqrt(A+B));

    // 			double Cij[bra::nc] __attribute__ ((aligned(16)));
    // 			double Cij_max = 0.0;

    // 			for(int j = 0, ij = 0; j < bra::B::nc; ++j) {
    // 			    for(int i = 0; i < bra::A::nc; ++i, ++ij) {
    // 				Cij[ij] = eAB*a(Ki,i)*b(Kj,j);
    // 				Cij_max = std::max(fabs(Cij[ij]), Cij_max);
    // 			    }
    // 			}
    // 			if (Cij_max*Ckl_max < cutoff) continue;

    // 			for(int kl = 0; kl < (c.nc*d.nc); ++kl) {
    // 			    for(int ij = 0; ij < bra::nc; ++ij) {
    // 				primitives.C[kl][ij + K*bra::nc] = Cij[ij]*Ckl[kl];
    // 			    }
    // 			}

    // 			Vector<3> rA = Vector<3>::center(ai, ri, aj, rj);
    // 			Vector<3> rAi = rA - ri;
    // 			Vector<3> rAB = rA - rB;

    // 			double rho = A*B/(A + B);
    // 			double X =  rho*rAB.inner();
    // 			Vector<N> W, t2;
    // 			rysq::roots<N>(X, t2, W);
    // 			t2 /= (A + B);

    // 			static const int LA = bra::A::L;
    // 			static const int LB = bra::B::L;

    // 			double *Ix = primitives.Ix(K);
    // 			double *Iy = primitives.Iy(K);
    // 			double *Iz = primitives.Iz(K); 

    // 			if (LB + d.L == 0) {
    // 			    rysq::recurrence<bra::L,N>(c.L + d.L, A, B, rAB, rAi, rBk, 
    // 						       t2, W, Ix, Iy, Iz);
    // 			}
    // 			else {
    // 			    F *Gx = primitives.template Gx<F>(K);
    // 			    F *Gy = primitives.template Gy<F>(K);
    // 			    F *Gz = primitives.template Gz<F>(K);
    // 			    F *tmp = primitives.template transfer<F>();

    // 			    rysq::recurrence<bra::L,N>(c.L + d.L, A, B, rAB, rAi, rBk, 
    // 						       t2, W, Gx, Gy, Gz);

    // 			    rysq::transfer<LA,LB,N>(c.L, d.L, rij[0], rkl[0], Gx, Ix, tmp);
    // 			    rysq::transfer<LA,LB,N>(c.L, d.L, rij[1], rkl[1], Gy, Iy, tmp);
    // 			    rysq::transfer<LA,LB,N>(c.L, d.L, rij[2], rkl[2], Gz, Iz, tmp);
    // 			}

    // 			++K;

    // 			// contract primitives
    // 			if (K == primitives.Kmax) {
    // 			    primitives.K = K;
    // 			    apply(c, d, primitives, scale, transform);
    // 			    K = 0;
    // 			}

    // 		    }
    // 		}
    // 	    }
    // 	}

    // 	// Final contraction
    // 	primitives.K = K;
    // 	if (primitives.K) apply(c, d, primitives, scale, transform);

    // }

private:

    template<typename F>
    static void apply(rysq::type c, rysq::type d,
		      const eri_::Primitives<F> &primitives,
		      double scale, Transform &transform);//  {

    // 	static const int ldN = (sizeof(F) > 8) ? N : N + N%2;
    // 	static const int Nij = ldN*(bra::A::L+1)*(bra::B::L+1);

    // 	const int Lc = abs(c);
    // 	const int Ld = abs(d);

    // 	const int dim2d = Nij*(Lc+1)*(Ld+1);

    // 	const int K = primitives.K;
    // 	const double* const *C = primitives.C;
    // 	const F *Ix = primitives.Ix(0);
    // 	const F *Iy = primitives.Iy(0);
    // 	const F *Iz = primitives.Iz(0);

    // 	int spk = (c < 0);
    // 	int spl = (d < 0);

    // 	const int c_first = shell(c).begin();
    // 	const int d_first = shell(d).begin();
    // 	const int c_last = shell(c).end() - 1;
    // 	const int d_last = shell(d).end() - 1;

    // 	for(int l = d_first, kl = 0; l <= d_last; ++l) {
    // 	    const int lx = (Lc+1)*LX[l];
    // 	    const int ly = (Lc+1)*LY[l];
    // 	    const int lz = (Lc+1)*LZ[l];

    // 	    const int lsp = (spl && l) << spk;

    // 	    for(int k = c_first; k <= c_last; ++k, ++kl) {
    // 		const double *Ckl = C[(spk && k) + lsp];

    // 		const int klx = Nij*(lx + LX[k]);
    // 		const int kly = Nij*(ly + LY[k]);
    // 		const int klz = Nij*(lz + LZ[k]);

    // 		static const int ni = bra::A::size;
    // 		static const int nj = bra::B::size;
    // 		double I[ni*nj] __attribute__((aligned(16))) = { 0.0 };
	    
    // 		int flags = 0;
    // 		double screen = 0.0;
    // 		quadrature::template apply<N>(flags, screen, K, Ckl, dim2d,
    // 					      &Ix[klx], &Iy[kly], &Iz[klz], 1.0, I);

    // 		transform(k - c_first, l - d_first, kl, I, scale);
    // 	    }
    // 	}

    // }

};


END_NAMESPACE(rysq, kernel)

#endif // _RYSQ_KERNEL_ERI_BRA_HPP_

