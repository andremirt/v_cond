#ifndef _RYSQ_KERNEL_ERI1_HPP_
#define _RYSQ_KERNEL_ERI1_HPP_

#include "kernel/quadrature.hpp"
#include <math.h>
#include <boost/mpl/int.hpp>
#include <rysq/core.hpp>
#include "vector.hpp"

#include "kernel/eri.hpp"
#include "cxx/utility/permute.hpp"

BEGIN_NAMESPACE(rysq, kernel)

template<typename T, size_t N>
struct align {
#define ALIGN_(n) ((A/sizeof(T) - (n)%(A/sizeof(T)))%(A/sizeof(T)))
    static const size_t A = 16;
    static const size_t value = ALIGN_(N);
    static size_t get(size_t M) { return ALIGN_(M); }
#undef ALIGN_
};

template<class bra_, int N>
struct Eri<bra_, boost::mpl::int_<N> > : public Eri<>
{
    typedef bra_ bra;
    typedef void ket;
    typedef kernel::Transform<bra> Transform;
    typedef Eri<>::Data Data;
    Eri(const Quartet<Shell> &quartet, Transform *transform)
	: quartet_(quartet),
	  transform_(transform)
    {
	primitives_.allocate<align>(quartet);
    }

    ~Eri() { delete transform_; }

    virtual void operator()(const Quartet<Center> &r, Data &data,
			    const Parameters &parameters) {
	double scale = rysq::SQRT_4PI5;
	double cutoff = parameters.cutoff/(scale*quartet_.K());
	quadrature::apply<bra,N,align>(quartet_, r[0], r[1], r[2], r[3],
				       scale, cutoff, primitives_,
				       (*transform_)(data));
    }

private:
    const Quartet<Shell> quartet_;
    Transform *transform_;
    quadrature::Primitives<double, double> primitives_;

};


END_NAMESPACE(rysq, kernel)

#endif // _RYSQ_KERNEL_ERI1_HPP_

