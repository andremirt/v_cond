#ifndef _RYSQ_ERI_IMPL_HPP_
#define _RYSQ_ERI_IMPL_HPP_

#include <boost/mpl/int.hpp>
#include <boost/utility/enable_if.hpp>

#include <rysq/eri.hpp>

#include "quartet.hpp"
#include "quadrature.h"
#include "normalize.h"
#include "util.h"
#include "kernels/functor.hpp"
#include "eri-eval.hpp"

#include "cxx/utility/permute.hpp"

namespace mpl = boost::mpl;
namespace meta = rysq::meta;


struct rysq::Eri::Impl {
    virtual ~Impl() {}
    virtual void operator()(const Quartet<Center> &centers,
			    double *I, const Parameters &parameters) = 0;
    template<type A, type B, type C, type D>
    static Impl* instance(const Quartet<Shell> &shells);
};

namespace rysq {
    namespace eri {


	template<class T, class Enable = void>
	struct Kernel;

	// template <type A_, type B_, type C_, type D_>
	// struct Kernel<meta::braket<A_,B_,C_,D_>,
	// 	      typename boost::disable_if<kernels::braket<A_,B_,C_,D_> >::type>
	// {
	//     static const bool value = false;
	// };

	template<class  braket, class enable = void>
	struct enable_if_specialized {};

	template <type A_, type B_, type C_, type D_>
	struct enable_if_specialized<meta::braket<A_,B_,C_,D_>,
				     typename 
				     boost::enable_if<kernels::braket<A_,B_,C_,D_>
						      >::type> {
	    typedef void type;
	};

	template <type A_, type B_, type C_, type D_>
	struct Kernel<meta::braket<A_,B_,C_,D_>, /*void*/
		      typename enable_if_specialized<
			  typename meta::sorted<A_,B_,C_,D_>::braket>:: type>
	: Eri::Impl {
	    static const bool value = true;
	    // typedef meta::braket<A_,B_,C_,D_>::braket braket;
	    typedef typename meta::sorted<A_,B_,C_,D_>::braket braket;
	    static const int transpose = meta::sorted<A_,B_,C_,D_>::value;

	    const Quartet<Shell> quartet;

	    Kernel(const Quartet<Shell> &quartet)
		: quartet(quartet.transposed(transpose)) { }

	    
	    static Impl* instance(const Quartet<Shell> &quartet) {
		return new Kernel(quartet);
	    }


	//     static void eval(const Quartet<Shell> &quartet,
	// 		     const Quartet<Center> &centers,
	// 		     const density_set &D6,
	// 		     const fock_set &F6,
	// 		     const Parameters &parameters) {
	// 	double I[braket::size] __attribute__ ((aligned(16))) = { 0.0 };
	// 	double scale = 1.0;
	// 	double cutoff = parameters.cutoff/(5e1*SQRT_4PI5);
	// 	eri::eval(quartet, centers, I, scale, cutoff);
	// 	Parameters p = parameters;
	// 	p.scale *= SQRT_4PI5;
	// 	fock::eval<braket>(D6, F6, I, p);
	//     }
	// };

	    void operator()(const Quartet<Center> &centers,
			    double *I, const Parameters &p) {
		double scale = SQRT_4PI5*p.scale;
		double cutoff = p.cutoff/(5e1*scale);
		eval(quartet, centers, I, scale, cutoff);
	    }

	    static void eval(const Quartet<Shell> &quartet,
			     const Quartet<Center> &centers,
			     double *I, double scale, double cutoff) {
		const double *ri, *rj, *rk, *rl;
		util::unpack(centers, ri, rj, rk, rl);
		cxx::utility::permute(ri, rj, rk, rl, quartet.shuffle_mask);
		eri::eval<braket>(quartet, ri, rj, rk, rl, I, scale, cutoff);
	    }

	};

	template<class braket, class Enable>
	struct Kernel {
	    static const bool value = false;
	    static Eri::Impl* instance(const Quartet<Shell> &quartet) {
		typedef typename braket::A A;
		typedef typename braket::B B;
		typedef typename braket::C C;
		typedef typename braket::D D;

		typedef typename braket::bra bra;
		typedef typename braket::ket ket;
		typedef mpl::int_<(braket::L)/2+1> N;

		typedef mpl::bool_<A::L == 0 && B::L != 0> t0;
		typedef mpl::bool_<C::L == 0 && D::L != 0> t1;
		typedef mpl::bool_<bra::L < ket::L && bra::L == 0> t01;

		// typedef mpl::bool_<A::value < B::value> t0;
		// typedef mpl::bool_<C::value < D::value> t1;
		// typedef mpl::bool_<((A::value == C::value && B::value < D::value) ||
		// 		    (A::value < C::value))> t01;

		typedef meta::state<B::type, A::type> BA_;
		typedef meta::state<A::type, B::type> AB_;
		typedef typename mpl::if_<t0, BA_, AB_>::type bra_;

		typedef meta::state<D::type, C::type> DC_;
		typedef meta::state<C::type, D::type> CD_;
		typedef typename mpl::if_<t1, DC_, CD_>::type ket_;

		typedef typename mpl::if_<t01, ket_, bra_>::type kBra;

		int transpose = Transpose::mask(t0::value, t1::value, t01::value);
		//if (braket::L < 5)std::cout << A::type << B::type << C::type  << D::type << "\n";
		return Kernel<kBra,N>::instance(quartet.transposed(transpose));
	    }
	};

    
	template<type A, type B, int N>
	struct Kernel<meta::state<A,B>, mpl::int_<N> >
	    : Eri::Impl, Quadrature::Transform {

	    typedef meta::state<A,B> bra;

	    const Quartet<Shell> quartet;
	    Quadrature::Primitives primitives;	    
	    Parameters parameters;
	    double scale, *eri;

	    Kernel(const Quartet<Shell> &quartet)
		: quartet(quartet),
		  primitives(this->quartet) {
	    }

	    static Eri::Impl* instance(const Quartet<Shell> &quartet) {
		return new Kernel(quartet);
	    }

	    void operator()(const Quartet<Center> &centers,
			    double *eri, const Parameters &parameters) {
		this->eri = eri;
		this->scale = parameters.scale*SQRT_4PI5;
		 double scale_ = 1;
		double cutoff_ = parameters.cutoff/(5e1*scale);
		Quadrature::eval<bra,N>(quartet, centers, scale_, cutoff_,
					primitives, *this);
	    }

	    struct Transform {
		static const int size = bra::size;
		double  scale, *eri;
		Transform(double scale,double *eri) : scale(scale), eri(eri) {}
		void operator()(int k, int l, int kl, const double *eri) {
		    // std::cout << kl << std::endl;
		    double *eri_ = this->eri + kl*size;
		    for (int i = 0; i < size; ++i) {
		        eri_[i] = scale*eri[i];
			// std::cout << eri[i] << std::endl;
		    }
		}
	    };

	    void operator()(const Quadrature::Primitives &primitives) {
		int K = primitives.K;
 		if (K == 0) return;
		const Shell &c = quartet[2];
		const Shell &d = quartet[3];
		Transform transform( scale,this->eri);
		// throw;
		Quadrature::apply<bra,N>(c, d, primitives, transform);
	    }
	    

	};

	// template <type A, type B, int N>
	// struct Kernel<meta::state<A,B>,  mpl::int_< N> >
	// : rysq::Eri::Impl, Quadrature::Transform {

	// };

    }
}


namespace rysq {

    template<type A, type B, type C, type D>
    inline Eri::Impl* Eri::Impl::instance(const Quartet<Shell> &quartet) {
	typedef meta::braket<A,B,C,D> braket;
	return eri::Kernel<braket>::instance(quartet);
    }

}


#endif /* _RYSQ_ERI_HPP_ */

