#include <boost/thread/tss.hpp>

#include "externals/cuda/host.hpp"
#include "externals/cuda/types.hpp"

#include <rysq/cuda.hpp>
#include "cuda/detail.hpp"
#include "index.hpp"

using namespace rysq::cuda;
namespace host = cuda::host;

namespace {

    using rysq::State;

    size_t coefficients(const State &bra, const State &ket, double2 *AB) {
	// size_t L = bra.L() + ket.L();
	// ket
#define COEFFICIENT(shell,index) ((shell)(index,(shell).is_hybrid()))
	for(int l = 0, ijkl = 0; l < ket[1].K; ++l) {
	    for(int k = 0; k < ket[0].K; ++k) {
		double B = ket[0](k) + ket[1](l);
		double Ckl = COEFFICIENT(ket[0], k)*COEFFICIENT(ket[1], l);
		// bra
		for(int j = 0; j < bra[1].K; ++j) {
		    for(int i = 0; i < bra[0].K; ++i, ++ijkl) {
			// compute values
			double A = bra[0](i) + bra[1](j);
			double Cij = COEFFICIENT(bra[0], i)*COEFFICIENT(bra[1], j);
			// double Cij = bra[0](i,0)*bra[1](j,0);
			double C = Cij*Ckl;
			double AB2 = 1.0/(sqrt(A+B)*A*B);
			AB[ijkl].x = C*AB2;
			AB[ijkl].y = 1.0/(A + B);
			// std::cout << AB[ijkl].x << AB[ijkl].y << "\n";
			// if (L == 0) AB[ijkl].y *= A*B;
			// if (L == 1) AB[ijkl].y *= A;
		    }
		}
	    }
	}
	return bra.K()*ket.K()*sizeof(double2);
#undef COEFFICIENT
    }

    size_t exponents(const State &state, double2 *dest) {
	for(int j = 0, ij = 0; j < state[1].K; ++j) {
	    for(int i = 0; i < state[0].K; ++i, ++ij) {
		dest[ij].x = state[0](i);
		dest[ij].y = state[1](j);
	    }
	}
	return state.K()*sizeof(double2);
    }

    size_t sp_coefficients(const rysq::Shell &shell, double *Cps) {
	if (!shell.is_hybrid()) return 0;
	for (int i = 0; i < shell.K; ++i) {
	    Cps[i] = shell(i,0)/shell(i,1);
	}
	return shell.K*sizeof(double);
    }

    size_t sp_coefficients(const State &state, double *Cps) {
	double *Cps0 = Cps;
	double *Cps1 = Cps + (state[0].type < 0)*state[0].K;
	for(int j = 0, ij = 0; j < state[1].K; ++j) {
	    double Cj = state[1](j,0);
	    for(int i = 0; i < state[0].K; ++i, ++ij) {
		double Ci = state[0](i,0);
		if (state[0].type < 0) Cps0[ij] = (state[0](i,1))/Ci;
		if (state[1].type < 0) Cps1[ij] = (state[1](j,1))/Cj;
	    }
	}
	return (state.K()*state.hybrid())*sizeof(double);
    }

    struct Index {
	static const Index& object() {
	    if (!Index::instance_.get()) {
		Index::instance_ .reset(new Index());
	    }
	    return *Index::instance_;
	}
	typedef host::array_ref<const ushort3> array_ref;
	Index() {
	    index host;
	    map_ = host.map();
	    const index::vector_type &data = host.data();
	    size_t size = data.size();
	    index_.reserve(size);
	    host::copy(data, index_);
	}
	array_ref find(rysq::Quartet<rysq::Shell> quartet) const {
	    map_type::const_iterator it = map_.find(index::key(quartet));
	    size_t size = quartet.size();
	    if (it == map_.end())
		return array_ref(host::wrap<ushort3>(NULL),  size);
	    // std::cout << index_.data().data() << " "
	    // 	  << (*it).second << " " << size << std::endl;
	    return array_ref(index_ + (*it).second, size);
	}

    private:
	static boost::thread_specific_ptr<Index> instance_;
	// static const Index *instance_;
	typedef rysq::Index<ushort3> index;
	typedef index::map_type map_type;
	map_type map_;
	host::array<ushort3> index_;
    };
    boost::thread_specific_ptr<Index> Index::instance_;

}

rysq::cuda::detail::Quartet::Quartet(const rysq::Quartet<Shell> &quartet)
    : base(quartet) {
    index2d_ = ::Index::object().find(quartet);

    const rysq::Quartet<Shell>::Bra &bra = quartet.bra();
    const rysq::Quartet<Shell>::Ket &ket = quartet.ket();
    size_t size = 0;

    size = (bra.K()*ket.K())*sizeof(double2);
    size += (bra.K() + ket.K())*sizeof(double2);
    for (int i = 0; i < 4; ++i)
	size += quartet[i].is_hybrid()*quartet[i].K*sizeof(double);

    char *data = new char[size];
    size_t bytes = 0;
    bytes += exponents(bra, (double2*)(data + bytes));
    bytes += exponents(ket, (double2*)(data + bytes));
    bytes += coefficients(bra, ket, (double2*)(data + bytes));

    for (int i = 0; i < 4; ++i)
	bytes += sp_coefficients(quartet[i], (double*)(data + bytes));

    assert(bytes == size);

    this->data_.reserve(size);
    host::copy(data, this->data_.data(), size);
    delete[] data;
}


