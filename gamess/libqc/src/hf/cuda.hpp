#ifndef _HF_CUDA_HPP_
#define _HF_CUDA_HPP_

#include <boost/thread.hpp>
#include <boost/multi_array.hpp>

#include "hf/hf.hpp"
#include "hf/thread.hpp"

#include "adapter/rysq.hpp"
#include <rysq.hpp>

#if RYSQ_CUDA
#define HF_CUDA 1

namespace hf {

    struct cuda {

	static std::vector<int> all_devices() {
	    return rysq::cuda::get_devices();
	}

	struct fock_thread_group {
	    boost::thread_group group_;
	    template<class B, class C, class S, class D, class F, class Q>
	    fock_thread_group(const B &blocks, const C &centers, const S &screen,
			      const D &density, lockable<F> &fock, Q &queue,
			      const std::vector<int> &devices = cuda::all_devices()) {
		for (size_t i = 0; i < devices.size(); ++i) {
		    group_.add_thread(fock_thread::create(blocks, centers, screen,
							  density, fock, queue,
							  devices[i]));
		}
	    }
	    void join_all() { group_.join_all(); }
	};

	struct fock_thread {

	    template<class B, class C, class S, class D, class F, class Q>
	    static boost::thread* create(const B &blocks, const C &centers,
					 const S &screen,
					 const D &density, lockable<F> &fock,
					 Q &queue, int device = 0) {
		using boost::ref;
		return new 
		    boost::thread(fock_thread::run<B, C, S, D, F, Q>,
				  blocks, ref(centers), ref(screen), 
				  ref(density), ref(fock), ref(queue), device);
	    }

	    template<class B, class C, class S, class D, class F, class Q>
	    fock_thread(const B &blocks, const C &centers, const S &screen,
			const D &density, lockable<F> &fock, Q &queue,
			int device = 0)
		: thread_(fock_thread::run<B, C, S, D, F, Q>,
			  blocks, boost::ref(centers), boost::ref(screen),
			  boost::ref(density), boost::ref(fock), boost::ref(queue),
			  device)
	    {
	    } 
	    void join() { thread_.join(); }

	private:
	    boost::thread thread_;
	    template<class B, class C, class S, class D, class F, class Q>
	    static void run(const B &blocks, const C &centers, const S &screen,
			    const D &density, lockable<F> &fock, Q &queue,
			    int device) {

		rysq::cuda::set_device(device);

		boost::array<size_t,2> shape = {{ density.array()[0].size(),
						  density.array()[1].size() }};

		density_array_type density_array(shape);
		fock_array_type fock_array(shape);

		for (size_t j = 0; j < shape[1]; ++j) {
		    for (size_t i = 0; i < shape[0]; ++i) {
			size_t size1 = density.m(i,j).size1();
			size_t size2 = density.m(i,j).size2();
			size_t block1 = density.m(i,j).block1();
			size_t block2 = density.m(i,j).block2();
			boost::array<int,2> ij = {{ i, j}};

			fock_array(ij).resize(size1, size2, block1, block2);
			density_array(ij).resize(size1, size2, block1, block2);

			fock_array(ij).clear();
			density_array(ij).assign(density.m(i,j));
		    }
		}

		rysq::cuda::Centers cuda_centers(centers);
		quartet_vector quartets;

		timer t;

		while (true) {
		    typename Q::value_type task;
		    try { task = queue.pop(); }
		    catch(std::exception&) { break; }

		    if (!run_task(blocks, cuda_centers, screen,
				  density_array, fock_array, task, quartets)) {
			queue.backlog().push(task);
		    }
		}

		std::cout << "GPU thread time: " << t << std::endl;

		fock.lock();
		for (size_t j = 0; j < fock_array[1].size(); ++j) {
		    for (size_t i = 0; i < fock_array[0].size(); ++i) {
			fock.data().m(i,j) += fock_array[i][j].copy();
		    }
		}
		fock.unlock();
		

	    }

	    typedef boost::multi_array<rysq::cuda::density_matrix,2> density_array_type;
	    typedef boost::multi_array<rysq::cuda::fock_matrix,2> fock_array_type;

	    struct quartet_vector {
		std::vector<rysq::Int4> host;
		rysq::cuda::Quartets cuda;
	    };


	    template<class B, class S, class T, class Q>
	    static bool run_task(const B blocks, const rysq::cuda::Centers &centers,
				 const S &screen,
				 density_array_type &density_array,
				 fock_array_type &fock_array,
				 const T &task, Q &quartets) {


		BOOST_AUTO(ib, &blocks[task[0]]);
		BOOST_AUTO(kb, &blocks[task[1]]);
		BOOST_AUTO(jb, &blocks[task[2]]);
		BOOST_AUTO(lb, &blocks[task[3]]);
		const BOOST_TYPEOF(ib) block[] = { ib, jb, kb, lb };

		if (!screen.test(block)) return true;

		typedef adapter::rysq::Shell Shell;
		Shell a(block[0]->shell());
		Shell b(block[1]->shell());
		Shell c(block[2]->shell());
		Shell d(block[3]->shell());

		rysq::cuda::Fock f(rysq::Quartet<rysq::Shell>(a,b,c,d));
		if (!f) return false;

		quartets.host.clear();
		screen(block, quartets.host);
		// if (quartets.host.size() == 0) return true;

		int base[4];
		for (int i = 0; i < 4; ++i)
		    base[i] = block[i]->firstShell();

		rysq::hf::matrix_ptr_set<rysq::cuda::fock_matrix> fock(base);
		rysq::hf::matrix_ptr_set<rysq::cuda::density_matrix> density(base);

		foreach (const size_t (&index)[2], rysq::index_list()) {
		    size_t i = index[0], j = index[1];
		    size_t im = std::distance(&blocks[0], block[i]);
		    size_t jm = std::distance(&blocks[0], block[j]);
		    // std::cout << im << " " << jm << std::endl;
		    fock.get(i,j) = &fock_array[im][jm];
		    density.get(i,j) = &density_array[im][jm];
		}

		quartets.cuda.assign(quartets.host);
		rysq::Parameters parameters(0, 1.0, 4.0, 5.0e-11);
		f(centers, quartets.cuda, density, fock, parameters);

		return true;
	    }
	    // }
	    // 	run_task(blocks, centers, screen, density, fock.data(), task, quartets);
	    // }
	    // while(true) {
	    // 	try { task = queue.backlog().pop(); }
	    // 	catch(...) { break; }
	    // 	run_task(blocks, centers, screen, density, fock.data(), task, quartets);
	    // }

	};
    };

    template<typename Screen>
    struct Cuda {
	typedef boost::multi_array<rysq::cuda::fock_matrix,2> fock_array_type;
	Cuda(const Basis &basis, const matrix::meta_matrix<BlockMatrix> &D,
	     const Screen &screen) :
	    num_blocks_(basis.blocks().size()),
	    centers_(basis.centers()),
	    density_array_(boost::extents[num_blocks_][num_blocks_]),
	    fock_array_(boost::extents[num_blocks_][num_blocks_]),
	    screen_(screen)
	{
	    for (size_t j = 0; j < num_blocks_; ++j) {
		for (size_t i = 0; i < num_blocks_; ++i) {
		    size_t size1 = D.m(i,j).size1();
		    size_t size2 = D.m(i,j).size2();
		    size_t block1 = D.m(i,j).block1();
		    size_t block2 = D.m(i,j).block2();
		    boost::array<int,2> ij = {{ i, j}};

		    fock_array_(ij).resize(size1, size2, block1, block2);
		    density_array_(ij).resize(size1, size2, block1, block2);

		    fock_array_(ij).clear();
		    density_array_(ij).assign(D.m(i,j));
		}
	    }
	}

	struct Thread : boost:: noncopyable{
	    const rysq::cuda::Centers &centers_;
	    adapter::rysq::Shell a_, b_, c_, d_;
	    rysq::cuda::Fock fock_;
	    template<class iterator>
	    Thread(const rysq::cuda::Centers &centers,
		   const iterator (&blocks)[4])
		: centers_(centers),
		  a_(blocks[0]->shell()),
		  b_(blocks[1]->shell()),
		  c_(blocks[2]->shell()),
		  d_(blocks[3]->shell()),
		  fock_(rysq::Quartet<rysq::Shell>(a_,b_,c_,d_))
	    {
		if (!fock_) throw std::exception();
	    }

	    void operator()(const std::vector<rysq::Int4> &quartets,
			    rysq::hf::matrix_ptr_set<rysq::cuda::density_matrix> D,
			    rysq::hf::matrix_ptr_set<rysq::cuda::fock_matrix> F) {
		rysq::cuda:: Quartets quartets_(quartets);
		rysq::Parameters parameters(0, 1.0, 4.0, 5.0e-11);
		fock_(centers_, quartets_, D, F, parameters);
	    }
	};


	// if (threads.first.timed_join(1)) {
	// 	delete threads.first;
	// 	threads.first = NULL;
	// 	std::swap(threads.first, threads.second);

	template<class iterator>
	bool thread(const iterator first, const iterator (&blocks)[4]) {
	    Thread *ct;
	    try { ct = new Thread(centers_, blocks); }
	    catch(...) { return false; }
	    // if (!ct->fock_) { delete ct;return false; }

	    std::vector<rysq::Int4> quartets;
	    quartets.clear();
	    screen_(blocks, quartets);
	    if (quartets.size() == 0) return true;

	    int base[4];
	    for (int i = 0; i < 4; ++i)
		base[i] = blocks[i]->firstShell();

	    rysq::hf::matrix_ptr_set<rysq::cuda::fock_matrix> F(base);
	    rysq::hf::matrix_ptr_set<rysq::cuda::density_matrix> D(base);

	    foreach (const size_t (&index)[2], rysq::index_list()) {
		size_t i = index[0], j = index[1];
		size_t im = std::distance(first, blocks[i]);
		size_t jm = std::distance(first, blocks[j]);
		// std::cout << im << " " << jm << std::endl;
		F.get(i,j) = &fock_array_[im][jm];
		D.get(i,j) = &density_array_[im][jm];
	    }

	    // boost::thread t(boost::ref(*ct),  quartets, D, F);
	    // t.join ();

	    (*ct)(quartets, D, F);
	    delete ct;

	    return true;
	}

	const fock_array_type& fock_array()const {return fock_array_; }
    private:
	size_t num_blocks_;
	rysq::cuda::Centers centers_;
	boost::multi_array<rysq::cuda::density_matrix,2> density_array_;
	fock_array_type fock_array_;
	const Screen &screen_;
    };

}

#endif /* RYSQ_CUDA */

#endif /* _HF_CUDA_HPP_ */
