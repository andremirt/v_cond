#ifndef _RYSQ_ERI_FACTORY_HPP_
#define _RYSQ_ERI_FACTORY_HPP_

#include <boost/preprocessor/seq/for_each_product.hpp>
#include <boost/preprocessor/seq/enum.hpp>
#include <boost/preprocessor/seq/elem.hpp>
#include "externals/cxx/namespace.hpp"

#include <rysq/core.hpp>
#include "kernel/eri.hpp"

BEGIN_NAMESPACE(rysq, eri)

struct Factory {
    template<class F>
    void new_(const Quartet<Shell> &quartet, F f) {
	type a = type(quartet[0]);
	type b = type(quartet[1]);
	type c = type(quartet[2]);
	type d = type(quartet[3]);

#define TYPES	(SP)(S)(P)(D)(F)

#define ERI(r, types) if (a == BOOST_PP_SEQ_ELEM(0, types) &&		\
			  b == BOOST_PP_SEQ_ELEM(1, types) &&		\
			  c == BOOST_PP_SEQ_ELEM(2, types) &&		\
			  d == BOOST_PP_SEQ_ELEM(3, types)) {		\
	    typedef kernel::eri<BOOST_PP_SEQ_ENUM(types)> kernel;	\
	    return new kernel(quartet, f.get<K::bra, K::ket>());	\
	} else								\

	BOOST_PP_SEQ_FOR_EACH_PRODUCT(ERI, (TYPES)(TYPES)(TYPES)(TYPES)) {
	    assert(0 && "unhandled exception");
	}
#undef ERI
#undef TYPES

    }

};

END_NAMESPACE(rysq, eri)

#endif /* _RYSQ_ERI_FACTORY_HPP_ */
