// #include <boost/preprocessor/seq/for_each_product.hpp>
// #include <boost/preprocessor/seq/to_tuple.hpp>

// #include <rysq/core.hpp>
// // #include "eri-impl.hpp"
// // #include "quadrature.h"
// #include "normalize.h"

// // #include "kernel/factory.hpp"

// using namespace rysq;

// Eri::Eri(const Quartet<Shell> &quartet) {

//     // kernel::new_<kernel::eri_::transform>(quartet);

//     this->size = quartet.size();

//     type a = type(quartet[0]);
//     type b = type(quartet[1]);
//     type c = type(quartet[2]);
//     type d = type(quartet[3]);

// #define _TYPES	(SP)(S)(P)(D)(F)

// #define _ERI(r, types) if (a == BOOST_PP_SEQ_ELEM(0, types) &&	\
// 			    b == BOOST_PP_SEQ_ELEM(1, types) &&	\
// 			    c == BOOST_PP_SEQ_ELEM(2, types) &&	\
// 			    d == BOOST_PP_SEQ_ELEM(3, types)) {	\
// 								\
// 	pImpl = Eri::Impl::instance<				\
// 	BOOST_PP_SEQ_ELEM(0, types),				\
// 	    BOOST_PP_SEQ_ELEM(1, types),			\
// 	    BOOST_PP_SEQ_ELEM(2, types),			\
// 	    BOOST_PP_SEQ_ELEM(3, types)>(quartet);		\
//     } else

//     BOOST_PP_SEQ_FOR_EACH_PRODUCT(_ERI, (_TYPES)(_TYPES)(_TYPES)(_TYPES)) {
// 	assert(0 && "unhandled exception");
//     }

// #undef _ERI
// #undef _TYPES

// }

// void Eri::operator()(const Quartet<Center> &centers,
// 		     double *I, const Parameters &parameters) {
//     (*pImpl)(centers, I, parameters);
// }

// void Eri::operator()(const std::vector<Center> &centers,
// 		     const std::vector<Int4> &quartets,
// 		      double *Eri, const Parameters &parameters) {

//     for (uint q = 0; q < quartets.size(); ++q) {

// 	int i, j, k, l;
// 	util::unpack(quartets[q], i, j, k, l);

// 	Quartet<Center> C(centers[i], centers[j], centers[k], centers[l]);
		
// 	Parameters p = parameters;
// 	//p.scale /= symmetry(i, j, k, l);
	
// 	//(*pImpl)(C, &Eri[q*size], p);
// 	(*pImpl)(C, &Eri[q* size], p);
//     }

// }


