#ifndef _RYSQ_CUDA_IMPL_HPP_
#define _RYSQ_CUDA_IMPL_HPP_

#include <cuda.h>
#include <cuda_runtime.h>

#include <boost/mpl/if.hpp>
#include <boost/mpl/bool.hpp>
// #include  <boost/type_traits/is_const.hpp>
#include <boost/array.hpp>
#include <boost/utility/enable_if.hpp>

#include <cmath>

#include <rysq/core.hpp>
#include <rysq/fock.hpp>
#include "rysq/cuda-memory.hpp"
#include "transpose.hpp"

#include "externals/cuda/host.hpp"

#include "externals/cxx/array.hpp"
#include "externals/cxx/utility/permute.hpp"
#include "externals/cxx/math/math.hpp"


namespace rysq {
    namespace cuda {

	namespace mpl = boost::mpl;

	namespace detail {

	    struct layout {
		template<typename T>
		T element(T i, T size1, T j, T size2);
	    };

	    template<typename T, bool mapped = false>
	    struct matrix_set {
		typedef rysq::index_list index_list;
		typedef detail::layout layout;


		typedef ushort size_array[4];

		size_array num_blocks, block_, base_index;
		// boost::array<ushort,4> num_blocks, block, base_index;
		T *data[6];

		template<class  T1, class T2, class T3>
		matrix_set(const T1 &num_blocks, const T2 &block, const T3 &index) {
		    for (int i = 0; i < 4; ++i)    {
		    	this->num_blocks[i] = num_blocks[i];
		    	this->block_[i] = block[i];
		    	this->base_index[i] = index[i];
		    }
		    
		}

		template<class S>
		explicit matrix_set(const S &s){
		    size_t size[] = { s.get(0,1)->size1(),
				      s.get(0,1)->size2(),
		    		      s.get(2,3)->size1(),
				      s.get(2,3)->size2() };
	    
		    size_t block[] = { s.get(0,1)->block1(),
				       s.get(0,1)->block2(),
		    		       s.get(2,3)->block1(),
				       s.get(2,3)->block2() };

		    for (int i = 0; i < 6; ++i)
		    	data[i] = s[i]->data();

		    for (int i = 0; i < 4; ++i)    {
		    	if (size[i]%block[i]) throw;
		    	this->num_blocks[i] = size[i]/block[i];
		    	this->block_[i] = block[i];
		    	this->base_index[i] = s.base_index[i];
		    }
		}

		const size_array& blocks()const {return block_; }

		__host__ __device__
		size_t data_size(size_t index) const {
		    return (cxx::math::multiply(num_blocks, size_t())*
			    tile_size(index));
		}

		__host__ __device__
		size_t  tile_size(size_t index) const {
		    return (block_[index_list::index1(index)]*
			    block_[index_list::index2(index)]);
		}

		template<class Quartet>
		__host__ __device__
		T* tile(size_t index, const Quartet &q){
		    size_t index1 =index_list::index1(index);
		    size_t index2 =index_list::index2(index);
		    return data[index] + block(index1, index2, q)*tile_size(index);
		}

		template<class Quartet>
		__host__ __device__
		const T* tile(size_t index, const Quartet &q) const {
		    return const_cast<matrix_set*>(this)->tile(index, q);
		}

		template<class Quartet>
		__host__ __device__
		int block(size_t index1, size_t index2, const Quartet &q) const {
		    return block(index1, index2, q, mpl::bool_<mapped>());
		    //return block(index1, index2, q, mpl::bool_<false>());
		}

		template<class Quartet>
		__host__ __device__
		int block(size_t index1, size_t index2, const Quartet &q,
			  const mpl::bool_<false>&) const {
		    return ((q[index1] - base_index[index1]) +
		    	    int(q[index2] - base_index[index2])*num_blocks[index1]);
		}

		template<class Quartet>
		__host__ __device__
		int block(size_t index1, size_t index2, const Quartet &q,
			  const mpl::bool_<true>&) const {
		    const size_t index3 = index_list::index3(index1, index2);
		    const size_t index4 = index_list::index4(index1, index2);
		    return (block(index1, index2, q, mpl::bool_<false>()) + 
			    block(index3, index4, q, mpl::bool_<false>())*
			    int(num_blocks[index1]*num_blocks[index2]));
		}
	    };

	    namespace host = ::cuda::host;

	    struct Quartet : rysq::Quartet<rysq::shell> {//: quartet_base {
		typedef rysq::Quartet<rysq::shell> base;
		Quartet() : base(0,0,0,0) {}
		Quartet(const rysq::Quartet<rysq::Shell> &quartet);
		operator rysq::Quartet<rysq::type>() const {
		    const Quartet &q = *this;
		    return rysq::Quartet<rysq::type>(q[0], q[1],q[2],q[3]);
		}
		size_t size() const {
		    const Quartet &q = *this;
		    return (q[0].size*q[1].size*q[2].size*q[3].size);
		}
		host::array_ref<const void> data() const { return data_; }
		host::array_ref<const ushort3> index2d() const { return index2d_; }
	    private:
		host::array<void> data_;
		host::array_ref<const ushort3> index2d_;
	    };

	    inline std::ostream& operator<<(std::ostream& stream, const Quartet &q) {
		return boost::operator<<(stream,q);
	    }


	    struct Quartets {
		typedef const Int4* const_pointer;
		const Int4 *data_;
		size_t size_;
		Quartets() {}
		Quartets(const Int4 *data, size_t size)
		    : data_(data), size_(size) {}
		const_pointer data()const {return data_; }
		size_t size() const { return size_; }
		__device__ Int4 get(size_t index, int permutation, void *shmem) const {
		    Int4 &tmp = *((Int4*)shmem);
		    tmp = data_[index];
		    // ::cuda::println(index4);
		    return  cxx::utility::permute(tmp, permutation);
		}
		__device__ Int4 operator[](size_t index) const {
		    return data_[index];
		}
	    };		



	    typedef const rysq::Center* Centers;

	    struct Eri {
		Eri(const rysq::Quartet<Shell> &quartet,
		    const rysq::Transpose &transpose);
		~Eri();
		void operator()(const Centers &centers,
					const Quartets &quartets,
					double *eri,
					const Parameters &parameters);
	    private: void *impl_;
	    };


	    struct Fock {
		typedef detail:: matrix_set<double, true>  mapped_fock_set;
		typedef detail:: matrix_set<double> density_set;
		typedef detail:: matrix_set<double>  fock_set;
		Fock(const rysq::Quartet<Shell> &quartet,
		     const rysq::Transpose &transpose);
		~Fock();
		void operator()(const detail::Centers &centers,
				const detail::Quartets &quartets,
				const density_set D, mapped_fock_set F,
				const Parameters &p);
	    private: void *impl_;
	    };

	}	    

    }
}

#endif /* _ RYSQ_CUDA_IMPL_HPP_ */
