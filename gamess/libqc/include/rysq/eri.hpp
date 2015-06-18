#ifndef _RYSQ_ERI_HPP_
#define _RYSQ_ERI_HPP_

#include "rysq/core.hpp"

namespace rysq {

    class Eri {
    public:
	class Impl;
	size_t size;
	Eri(const Quartet<Shell> &shells);
	void operator()(const Quartet<Center> &centers,
			double *I, const Parameters &parameters);
	void operator()(const std::vector<Center> &centers,
			const std::vector<Int4> &quartets,
			double *Eri,
			const Parameters &parameters = rysq::Parameters());
    private:
	Impl *pImpl;
    };

}


#endif /* _RYSQ_ERI_HPP_ */
