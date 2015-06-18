#include <boost/thread/tss.hpp>
#include <boost/typeof/typeof.hpp>
#define AUTO BOOST_AUTO

#include "rysq/cuda.hpp"
#include "cuda/detail.hpp"
#include "cuda/quartet.hpp"
#include "transpose.hpp"

#include "externals/cxx/foreach.hpp"
#include "externals/cxx/sugar/unpack.hpp"
#include "externals/cxx/math/math.hpp"

#include "externals/cuda/host.hpp"

#include "cuda/matrix.hpp"

using namespace rysq;
namespace host = ::cuda::host;

rysq::cuda::Fock::Fock(const rysq::Quartet<rysq::Shell> &quartet) {
    bool t0 = (!quartet[0].L && quartet[1].L);
    bool t1 = (!quartet[2].L && quartet[3].L);
    bool t01 = (!quartet.bra().L() && quartet.ket().L());
    Transpose transpose(t0, t1, t01);

    try { this->kernel_ = new cuda::detail::Fock(quartet, transpose); }
    catch (std::exception&) { this->kernel_ = NULL; }

    for (int i = 0; i < 4; ++i)	block_[i] = quartet[i].size;
}

rysq::cuda::Fock::~Fock() {
    delete this->kernel_;
}

void rysq::cuda::Fock::operator()(const Centers &centers,
				  const Quartets &quartets,
				  density_matrix_set D, fock_matrix_set F,
				  const Parameters &parameters) {

    boost::array<size_t,4> mapped_dims;
    boost::array<int,4> mapped_index;

    for (int i = 0; i < 4; ++i) {
	mapped_dims[i] = quartets.max[i] - quartets.min[i] + 1;
	mapped_index[i] = quartets.min[i] - F.base_index[i];
    }

    detail::matrix_set<double, true> mapped_fock(mapped_dims, block_, quartets.min);    
    size_t mapped_size = 0;
    for (int i = 0; i < 6; ++i) {
	mapped_size += mapped_fock.data_size(i);
    }

    static boost::thread_specific_ptr<rysq::cuda::array<double> > thread_mapped_data;
    if (!thread_mapped_data.get()) {
	thread_mapped_data.reset(new rysq::cuda::array<double>());
    }

    // std::cout << mapped_size* sizeof(double) << std::endl;
    BOOST_AUTO(&mapped_data, *thread_mapped_data);
    try {
	mapped_data.resize(mapped_size);
    }
    catch (std::exception&) {
	// 128 mb is not available
	if (mapped_size*sizeof(double) < 128*(1024*1024)) throw;
	// recursively half memory
	size_t half = quartets.size()/2;
	this->operator()(centers, Quartets(quartets.begin(), quartets.begin() + half),
			 D, F, parameters);
	this->operator()(centers, Quartets(quartets.begin() + half, quartets.end()),
	 		 D, F, parameters);
	return;
    }

    cuda::memset(mapped_data, char(0), mapped_data.size()* sizeof(double));

    for (size_t i = 0, offset = 0; i < 6; ++i) {
	mapped_fock.data[i] = mapped_data.data() + offset;
	offset += mapped_fock.data_size(i);
    }

    detail::matrix_set<double> density(D);
    detail::matrix_set<double> fock(F);
    

    (*kernel_)(detail::Centers(centers.data()),
    	       detail::Quartets(quartets.data(), quartets.size()),
    		density, mapped_fock, parameters);

    
    foreach (const size_t (&index)[2], index_list()){
    	size_t i = index[0], j = index[1];
    	size_t k = index_list::index3(i,j);
    	size_t l = index_list:: index4(i,j);
	double scale = (i + j == 1 || i + j == 5) ? 4 : -1;

    	size_t block_size = block_[i]*block_[j];
    	size_t size1 =  mapped_dims[i]*block_size, size2 = mapped_dims[j];
    	size_t ld = fock.num_blocks[i]*block_size;
    	double *data = F.get(i,j)->block(mapped_index[i],  mapped_index[j]);

    	rysq::matrix<double> fock_(size1, size2, data, ld);
    	rysq::matrix_data_array<double> 
    	    mapped_fock_(matrix_layout(size1, size2),
    			 mapped_fock.data[index_list::find(index)],
    			 mapped_dims[k]*mapped_dims[l]);
    	cuda::reduce(mapped_fock_, scale, fock_);
    }

}

