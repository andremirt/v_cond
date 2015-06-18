// /**
//  * @file   
//  * @brief  Fock JK matrix operator implementation
//  */

// #include <algorithm>

// #include "hf/work_queue.hpp"
// #include "hf/hf.hpp"
// #include "hf/thread.hpp"
// #include "runtime.hpp"

// #include "util/numeric.hpp"
// #include "matrix/meta.hpp"

// #include "hf/cuda.hpp"

// #include <boost/array.hpp>
// #include <boost/typeof/typeof.hpp>
// #include <boost/utility.hpp>
// #include <boost/math/special_functions/pow.hpp>
// #include <boost/date_time/posix_time/posix_time.hpp>
// #include <boost/lambda/lambda.hpp>

// using boost::math::pow;
// namespace lambda = boost::lambda;

// #define TYPEOF BOOST_TYPEOF
// #define AUTO BOOST_AUTO

// #ifdef HAVE_GA
// #include "distributed.hpp"
// #endif


// using namespace hf;

// // see hf.h
// void hf::fock(const Basis &basis, 
// 	      const matrix::meta_matrix<BlockMatrix> &D,
// 	      matrix::meta_matrix<BlockMatrix> &F,
// 	      const matrix::meta_matrix<Matrix> *Kmax,
// 	      const matrix::meta_matrix<Matrix> *Dmax,
// 	      double cutoff) {

//     timer t;

//     Screen<matrix::meta_matrix<Matrix> > screen(*Kmax, *Dmax, cutoff);

//     typedef rysq::block_matrix_adapter<double> matrix_adapter;
//     rysq::initialize();

//     work_queue queue(basis.blocks().size());

    

//     lockable<TYPEOF(F)> lockable_fock(F);

//     size_t num_threads = runtime::num_threads(fock_thread_group::max_threads());
//     fock_thread_group host(basis.blocks().begin(), basis.centers(), screen,
// 			   D, lockable_fock, queue, num_threads);

// #if HF_CUDA
//     {
// 	std::vector<int> devices = runtime::cuda_devices(cuda::all_devices());
// 	cuda::fock_thread_group cuda(basis.blocks().begin(), basis.centers(), screen,
// 				     D, lockable_fock, queue, devices);
// 	cuda.join_all();
//     }
// #endif

//     host.join_all();

//     if (!queue.backlog().empty()) {
//     	fock_thread host(basis.blocks().begin(), basis.centers(), screen,
//     			 D, lockable_fock, queue);
//     	host.join();
//     }

//     matrix::symmeterize(F); // symmetrize fock matrix

//     // surface.render(K);
//     //matrix::plot::surface(F);

//     return;
// }

