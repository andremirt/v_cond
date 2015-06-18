#ifndef _RYSQ_ERI_NEW_HPP_
#define _RYSQ_ERI_NEW_HPP_

#include <boost/preprocessor/seq/for_each_product.hpp>
#include <boost/preprocessor/seq/enum.hpp>
#include <boost/preprocessor/seq/elem.hpp>
#include "externals/cxx/namespace.hpp"

#include <rysq/core.hpp>
#include "kernel/eri.hpp"

BEGIN_NAMESPACE(rysq, kernel)




template<template<class, class> class T>
kernel::Eri<>* new_(const Quartet<Shell> &quartet) {
    type a = type(quartet[0]);
    type b = type(quartet[1]);
    type c = type(quartet[2]);
    type d = type(quartet[3]);

#define TYPES	(rysq::SP)(rysq::S)(rysq::P)(rysq::D)(rysq::F)

#define ERI(r, types) if (a == BOOST_PP_SEQ_ELEM(0, types) &&		\
			  b == BOOST_PP_SEQ_ELEM(1, types) &&		\
			  c == BOOST_PP_SEQ_ELEM(2, types) &&		\
			  d == BOOST_PP_SEQ_ELEM(3, types)) {		\
									\
	typedef typename						\
	    kernel::find<BOOST_PP_SEQ_ENUM(types)>::type kernel;	\
	return new kernel(quartet, new T<kernel::bra, kernel::ket>());	\
    } else								\

    BOOST_PP_SEQ_FOR_EACH_PRODUCT(ERI, (TYPES)(TYPES)(TYPES)(TYPES)) {
	throw std::runtime_error("invalid type");
    }
#undef ERI
#undef TYPES

}


END_NAMESPACE(rysq, eri)

#endif /* _RYSQ_ERI_NEW_HPP_ */
