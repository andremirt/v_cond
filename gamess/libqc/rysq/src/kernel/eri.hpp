#ifndef _RYSQ_KERNEL_ERI_HPP_
#define _RYSQ_KERNEL_ERI_HPP_


#include <boost/mpl/int.hpp>
#include <boost/mpl/if.hpp>
#include <boost/utility/enable_if.hpp>

#include "rysq/core.hpp"
#include "meta.hpp"
#include "externals/cxx/namespace.hpp"


BEGIN_NAMESPACE(rysq, kernel)

template<class bra, class ket = void, class enabled = void>
struct eri {
    static const bool value = false;
};

template<class bra = void, class ket = void>
struct Transform;

template<>
struct Transform<> {
    struct Data {};
};

template<class bra>
struct Transform<bra> {
    typedef typename Transform<>::Data Data;
    virtual ~Transform() {}
    virtual Transform& operator()(Data &data) = 0;
    virtual void operator()(int k,int l,int kl,
			    const double *Q, double scale) = 0;
};

template<class bra, class ket>
struct Transform {
    typedef typename Transform<>::Data Data;
    virtual ~Transform() {}
    virtual Transform& operator()(Data &data) = 0;
    virtual void operator()(const double *Q, double scale) = 0;
};

template<class A = void, class B = void, class enable = void>
struct Eri {
    static const bool value = false;
};

template<>
struct Eri<> {
    typedef Transform<>::Data Data;
    virtual ~Eri() {}
    virtual void operator()(const Quartet<Center> &r, Data &data,
			    const Parameters &parameters) = 0;
};


END_NAMESPACE(rysq, kernel)

#include "kernel/eri1.hpp"
#include "kernel/eri2.hpp"

BEGIN_NAMESPACE(rysq, kernel)

template<type T0, type T1, type T2, type T3>
struct find {
    typedef typename meta::state<T0,T1> bra;
    typedef typename meta::state<T2,T3> ket;
    //static const size_t L = bra::L + ket::L;
    typedef boost::mpl::int_<(bra::L + ket::L)/2 + 1> roots;
    typedef typename boost::mpl::if_<Eri<bra,ket>,
				     Eri<bra,ket>,
				     Eri<bra,roots>
				     >::type type;
};

END_NAMESPACE(rysq, kernel)

#endif /* _RYSQ_KERNEL_ERI_HPP_ */
