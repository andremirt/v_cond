#ifndef _HF_SCREEN_HPP_
#define _HF_SCREEN_HPP_

#include "hf/hf.hpp"

namespace hf {

    template<class M1, class M2>
    void Dmax(const Basis &basis, const M1 &D, M2 &Dmax) {
	for (int b = 0; b < basis.numShells(); ++b) {
	    int nj = basis(b).size();
	    int j0 = basis.firstShellFunction(b);

	    for (int a = 0; a <= b; ++a) {
		int ni = basis(a).size ();
		int i0 = basis.firstShellFunction(a);

		double q = 0.0;
		for (int j = j0; j < j0+nj; ++j) {
		    for (int i = i0; i < i0+ni; ++i) {
			q = std::max(q, fabs(D(i,j)));
		    }
		}
		Dmax(a,b) = q;
		Dmax(b,a) = q;

	    }
	}
    }

    template<class M1, class M2 = M1>
    class Screen {
    public:
	Screen(const M1 &Kmax, const M2 &Dmax, double cutoff) :
	    Kmax_(Kmax), Dmax_(Dmax), cutoff_(cutoff) {}

	template<typename  iterator, class C>
	void operator()(const iterator (&block)[4], C &quartets) const {
	    (*this)(*block[0], *block[1], *block[2], *block[3], quartets);
	}
	template<typename  iterator>
	bool test(const iterator (&block)[4]) const {
	    return test(*block[0], *block[1], *block[2], *block[3]);
	}

	bool test(const Basis::Block &a, const Basis::Block &b,
			const Basis::Block &c, const Basis::Block &d)  const {
	    std::vector<boost::array<bool,4> > quartets;
	    return apply<true>(a,b,c,d,  quartets);
	}
	template<class C>
	bool operator()(const Basis::Block &a, const Basis::Block &b,
			const Basis::Block &c, const Basis::Block &d,
			C &quartets)  const {
	    return apply<false>(a,b,c,d,  quartets);
	}

	template<bool test_only, class C>
	bool apply(const Basis::Block &a, const Basis::Block &b,
		   const Basis::Block &c, const Basis::Block &d,
		   C &quartets)  const {

	    // build vector of quartets
	    for (int l = d.firstShell(); l <= d.lastShell(); ++l) {
		for (int j = b.firstShell(); j <= b.lastShell(); ++j) {

		    int kfirst = (c == d) ? l : c.firstShell();
		    for (int k = kfirst; k <= c.lastShell(); ++k) {

			int ifirst = (a == b) ? j : a.firstShell();
			int ilast = (a == c) ? k + (j <= l) : a.lastShell() + 1;
			for (int i = ifirst; i < ilast; ++i) {
			    if (!test(i, j, k, l)) continue;
			    if (test_only) return true;
			    typename C::value_type quartet = {{i, j, k, l}};
			    quartets.push_back(quartet);
			}
		    }

		}
	    }
	    return false;
	}
	bool test(int i,int j,int k,int l)const {
	    return (max(Dmax_,i, j, k, l)*Kmax_(i,j)*Kmax_(k,l) > cutoff_);
	}
	template<typename T>
	static T max(const T &a, const T &b, const T &c, const T &d) {
	    return std::max(std::max(a,b), std::max(c,d));
	}
	template<class M>
	static double max(const M &Dmax, int i, int j, int k, int l) {
	    return std::max((4.0*std::max(Dmax(i,j), Dmax(k,l))),
			    (max(Dmax(i,k), Dmax(i,l), Dmax(j,k), Dmax(j,l))));
	}
    private:
	const M1 &Kmax_;
	const M2 &Dmax_;
	double cutoff_;
    };

}

#endif /* _HF_SCREEN_HPP_ */
