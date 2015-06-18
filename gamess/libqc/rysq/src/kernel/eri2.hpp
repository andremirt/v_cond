#ifndef _RYSQ_KERNEL_ERI2_HPP_
#define _RYSQ_KERNEL_ERI2_HPP_

#include "kernel/eri.hpp"
#include "kernel/quadrature.hpp"
#include "kernel/quadrature2-impl.hpp"
#include "transpose.hpp"

BEGIN_NAMESPACE(rysq, kernel)

template<type T0, type T1, type T2, type T3>
struct Eri<meta::state<T0, T1>, meta::state<T2, T3>,
	   typename quadrature::impl<
	       typename meta::sorted<T0,T1,T2,T3>::braket>::enable >
    : public Eri<>
{
    static const bool value = true;
    typedef meta::state<T0, T1> bra;
    typedef meta::state<T2, T3> ket;
    typedef typename meta::sorted<T0,T1,T2,T3>::braket braket;
    static const int transpose_value = meta::sorted<T0,T1,T2,T3>::value;

    typedef kernel::Transform<bra,ket> Transform;
    typedef Eri<>::Data Data;

    Eri(const Quartet<Shell> &quartet, Transform *transform)
	: quartet_(transpose(quartet, transpose_value)),
	  transform_(transform) {
	// std::cout << quartet << transpose(quartet, transpose_value) << std::endl;	
    }

    ~Eri() { delete transform_; }

    void operator()(const Quartet<Center> &r,
		    Data &data,
		    const Parameters &parameters) {
	const Quartet<Center> r_ = transpose(r, transpose_value);
	double scale = rysq::SQRT_4PI5;
	double cutoff = parameters.cutoff/(scale*quartet_.K());
	double Q[braket::size]  __attribute__ ((aligned(16))) = { 0.0 };
	// std::cout << quartet_ << transpose_value<< std::endl;	
	quadrature::apply<braket>(quartet_, r_[0], r_[1], r_[2], r_[3],
				  Q, transpose_value, cutoff);
	((*transform_)(data))(Q, scale);
    }

private:
    const Quartet<Shell> quartet_;
    Transform *transform_;
};

END_NAMESPACE(rysq, kernel)

#endif /* _RYSQ_KERNEL_ERI2_HPP_ */
