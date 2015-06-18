#ifndef _RYSQ_CUDA_HPP_
#define _RYSQ_CUDA_HPP_

#include <vector>

// #include <boost/numeric/ublas/matrix.hpp>
#include "rysq/core.hpp"
#include "rysq/eri.hpp"
#include "rysq/fock.hpp"


#include "rysq/cuda-memory.hpp"
#include "rysq/cuda-matrix.hpp"

namespace rysq {
    namespace cuda {

	std::vector<int> get_devices();
	void set_device(size_t device = 0);
	void thread_exit();

	struct thread {
	    thread(size_t device = 0) { set_device(device); }
	    ~thread() { thread_exit(); }
	};

	namespace detail {
	    class Eri;
	    class Fock;
	}

	typedef cuda::array<Center> Centers;

	struct Quartets : cuda::array<Int4> {
	    Quartets() {
		min.assign(1);
		max.assign(0);
	    }
	    template<typename T>
	    explicit Quartets(T begin, T end)
		: vector_(begin, end) {
		initialize();
	    }
	    explicit Quartets(const std::vector<Int4> &centers)
		: vector_(centers) {
		initialize();
	    }
	    void assign(const std::vector<Int4> &centers) {
		vector_ = centers;
		initialize();
	    }
	    typedef std::vector<Int4>::const_iterator const_iterator;
	    const_iterator begin() const { return vector_.begin(); }
	    const_iterator end() const { return vector_.end(); }
	    boost::array<int,4> max, min;
	private:
	    typedef cuda::array<Int4> base_type;
	    std::vector<Int4> vector_;
	    void initialize() {
		max.assign(0);
		min.assign(1);
		if (!vector_ . empty()) {
		    max = vector_.front();
		    min = max;
		}
		typedef std::vector<Int4>::const_iterator iterator;
		for (iterator it = vector_.begin(); it < vector_.end(); ++ it) {
		    // std:: cout << *it << std::endl;
		    for (int i = 0; i < 4; ++i) {
			max[i] = std::max(max[i], (*it)[i]);
			min[i] = std::min(min[i], (*it)[i]);
		    }
		}
		base_type::assign(vector_);
	    }
	};

	struct Eri {
	    Eri(const rysq::Quartet<rysq::Shell> &shells);
	    ~Eri();
	    operator bool() const { return kernel_ != NULL; }
	    void operator()(const cuda::Centers &centers,
			    const cuda::Quartets &quartets,
			    cuda::array <double> &Eri,
			    const Parameters &parameters);
	private:
	    detail::Eri *kernel_;
	};

	typedef cuda::block_matrix<double> density_matrix;
	typedef cuda::block_matrix<double> fock_matrix;

	struct Fock {
	    typedef hf::matrix_ptr_set < fock_matrix> fock_matrix_set;
	    typedef hf::matrix_ptr_set < density_matrix> density_matrix_set;
	    Fock(const rysq::Quartet<rysq::Shell> &shells);
	    ~Fock();
	    operator bool() const { return kernel_ != NULL; }
	    // operator bool() const { return  kernel_ != NULL; }
	    void operator()(const cuda::Centers &centers,
			    const cuda::Quartets &quartets,
			    density_matrix_set D, fock_matrix_set F,
			    const Parameters &parameters);
	private:
	    detail::Fock *kernel_;
	    size_t block_[4];
	};
	    

    }
}

#endif /* _RYSQ_CUDA_HPP_ */
