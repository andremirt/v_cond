#ifndef _RYSQ_FOCK_IMPL_HPP_
#define _RYSQ_FOCK_IMPL_HPP_

#include <boost/mpl/int.hpp>
#include <boost/mpl/if.hpp>
#include <boost/utility/enable_if.hpp>

#include "externals/cxx/namespace.hpp"


#include <rysq/fock.hpp>

#include "fock-matrix.hpp"
#include "quartet.hpp"
#include "eri-impl.hpp"
#include "meta.hpp"
#include "transpose.hpp"
#include "fock-eval.hpp"


struct rysq::Fock::Impl {
    virtual ~Impl() {};
    virtual void operator()(const Quartet<Center> &centers,
			    const density_set &D,
			    const fock_set &F,
			    const Parameters &parameters) = 0;

    template<type A, type B, type C, type D>
    static Impl* instance(const Quartet<Shell> &shells);

};

BEGIN_NAMESPACE(rysq, fock)

namespace mpl = boost::mpl;

template<class braket, class Enable = void>
struct Kernel;

template<type A_, type B_, type C_, type D_>
struct Kernel<meta::braket<A_,B_,C_,D_>, /*void*/
	      typename boost::enable_if<
		  eri::Kernel<typename meta::sorted<A_,B_,C_,D_>::braket>
		  >::type> : Fock::Impl {

    typedef meta::braket<A_,B_,C_,D_> braket;
    typedef eri::Kernel<typename meta::sorted<A_,B_,C_,D_>::braket> eri;
    static const int transpose = meta::sorted<A_,B_,C_,D_>::value;

    const Quartet<Shell> quartet;

    Kernel(const Quartet<Shell> &quartet)
	: quartet(quartet.transposed(transpose)) { }

    static Fock::Impl* instance(const Quartet<Shell> &quartet) {
	return new Kernel(quartet);
    }

    void operator()(const Quartet<Center> &centers,
		    const density_set &D6,
		    const fock_set &F6,
		    const Parameters &parameters) {
	eval(quartet, centers, D6, F6, parameters);
    }

    static void eval(const Quartet<Shell> &quartet,
		     const Quartet<Center> &centers,
		     const density_set &D6,
		     const fock_set &F6,
		     const Parameters &parameters) {
	double I[braket::size] __attribute__ ((aligned(16))) = { 0.0 };
	double scale = 1.0;
	double cutoff = parameters.cutoff/(5e1*SQRT_4PI5);
	eri::eval(quartet, centers, I, scale, cutoff);
	Parameters p = parameters;
	p.scale *= SQRT_4PI5;
	fock::eval<braket>(D6, F6, I, p);
    }
};


template<class braket, class Enable>
struct Kernel {
    static Fock::Impl* instance(const Quartet<Shell> &quartet) {
	typedef typename braket::A A;
	typedef typename braket::B B;
	typedef typename braket::C C;
	typedef typename braket::D D;

	typedef typename braket::bra bra;
	typedef typename braket::ket ket;
	typedef mpl::int_<(braket::L)/2+1> N;

	// 		typedef mpl::bool_<A::L == 0 && B::L != 0> t0;
	// 		typedef mpl::bool_<C::L == 0 && D::L != 0> t1;
	// 		typedef mpl::bool_<bra::L < ket::L && bra::L == 0> t01;

	typedef mpl::bool_<((A::value == C::value && B::value < D::value) ||
			    (A::value < C::value))> t01;

	typedef meta::transpose<A, B, (A::value < B::value)> bra_;
	typedef meta::transpose<C, D, (C::value < D::value)> ket_;

	typedef typename mpl::if_<t01,
	    typename ket_::type, typename bra_::type>::type kBra;

	int transpose = Transpose::mask(bra_::value, ket_::value, t01::value);
	return NULL;
    }
};


    
template<type A, type B, int N>
struct Kernel<meta::state<A,B>, mpl::int_<N> >
    : Fock::Impl, Quadrature::Transform {

    typedef meta::state<A,B> bra;

    const Quartet<Shell> quartet;
    Quadrature::Primitives primitives;	    
    Parameters parameters;
    fock::Matrix matrix;

    Kernel(const Quartet<Shell> &quartet)
	: quartet(quartet),
	  primitives(this->quartet),
	  matrix(this->quartet) {
    }

    static Fock::Impl* instance(const Quartet<Shell> &quartet) {
	return new Kernel(quartet);
    }

    void operator()(const Quartet<Center> &centers,
		    const density_set &D,
		    const fock_set &F,
		    const Parameters &parameters) {
	matrix.set(D);

	double scale = 1.0;
	double cutoff = parameters.cutoff/(5e1*SQRT_4PI5);
	Quadrature::eval<bra,N>(quartet, centers, scale, cutoff, primitives, *this);

	scale = parameters.scale*SQRT_4PI5;
	double scale2 = parameters.scale2*scale;
	matrix.get(scale, scale2, F);
    }

    struct Transform {
	static const int ni = bra::A::size;
	static const int nj = bra::B::size;
	density_set D;
	fock_set F;
	Transform(const density_set &D, const fock_set &F)
	    : D(D), F(F) {}
	void operator()(int k, int l, int kl, const double *I) {
	    fock::eval<ni,nj>(k, l, kl, D, F, I);
	}
    };

    void operator()(const Quadrature::Primitives &primitives) {
	int K = primitives.K;
	if (K == 0) return;
	const Shell &c = quartet[2];
	const Shell &d = quartet[3];
	Transform transform(matrix.D(), matrix.F());
	Quadrature::apply<bra,N>(c, d, primitives, transform);
    }


};

END_NAMESPACE(rysq,fock)

namespace rysq {

    typedef Fock::density_set density_set;
    typedef Fock::fock_set fock_set;

    template<class bra, class ket = void>
    struct Transform {
	static const size_t size = bra::size*ket::size;
	density_set D;
	fock_set F;
	Parameters parameters;
	Transform(const density_set &D, const fock_set &F,
		  const Parameters &parameters)
	    : D(D), F(F), parameters(parameters) {}
	void operator()(double scale, const double (&Q)[size], size_t num = size) {
	    Parameters p = parameters;
	    p.scale *= SQRT_4PI5;
	    fock::eval<bra,ket>(D, F, Q, p);
	}
    };

    template<class bra>
    struct Transform<bra,void> {
	static const int ni = bra::A::size;
	static const int nj = bra::B::size;
	density_set D;
	fock_set F;
	Parameters parameters;
	Transform(const density_set &D, const fock_set &F,
		  const Parameters &parameters)
	    : D(D), F(F), parameters(parameters) {}
	void operator()(int k, int l, int kl, const double *I) {
	    fock::eval<ni,nj>(k, l, kl, D, F, I);
	}
	void operator()(const double Q[bra::size]) {
	    double scale = 1.0;
	    Parameters p = parameters;
	    p.scale *= SQRT_4PI5;
	    fock::eval<bra>(D, F, Q, p);
	}
    };
	    



    template<type A, type B, type C, type D>
    inline Fock::Impl* Fock::Impl::instance(const Quartet<Shell> &quartet) {
	typedef meta::braket<A,B,C,D> braket;
	return fock::Kernel<braket>::instance(quartet);
    }


}

#endif /* _RYSQ_FOCK_HPP_ */
