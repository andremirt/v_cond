#ifndef _HF_THREAD_HPP_
#define _HF_THREAD_HPP_

#include <vector>
#include <boost/array.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread.hpp>

#include <boost/typeof/typeof.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

// typedef boost::posix_time::ptime ptime;
// typedef boost::posix_time::time_duration time_duration;
// // typedef boost::posix_time::microsec_clock microsec_clock;

#include "hf/hf.hpp"

#include <rysq.hpp>
#include "adapter/rysq.hpp"

namespace hf {


    struct timer {
	boost::posix_time::ptime start_;
	typedef boost::posix_time::microsec_clock microsec_clock;
	timer() : start_(microsec_clock::universal_time()) {}
	boost::posix_time::time_duration duration() const {
	    return microsec_clock::universal_time() - start_;
	}
	long total_seconds() const {
	    return duration().total_seconds();
	}
    };
    std::ostream& operator << (std::ostream &os, const timer &t) {
	return os << boost::posix_time::to_simple_string(t.duration());
    }

    template<typename T>
    class lockable {
    public:
	typedef T& reference;
	typedef const T& const_reference;
	explicit lockable(T &data) : data_(data) {}
	void lock() const { mutex_.lock(); }
	void unlock() const { mutex_.unlock(); }
	reference data(){return data_; }
	const_reference data()const {return data_; }
    private:
	mutable boost:: mutex mutex_;
	T &data_;
    };

    template<typename T>
    lockable<T> make_lockable(T &data) {
	return lockable <T>(data);
    }

    class fock_thread {
    public:

	template<class B, class C, class S, class D, class F, class Q>
	static boost::thread* create(const B &blocks, const C &centers,
				     const S &screen,
				     const D &density, lockable<F> &fock,
				     Q &queue) {
	    using boost::ref;
	    return new 
		boost::thread(fock_thread::run<B, C, S, D, F, Q>,
			      blocks, ref(centers), ref(screen), 
			      ref(density), ref(fock), ref(queue));
	}

	boost::thread thread_;
	template<class B, class C, class S, class D, class F, class Q>
	fock_thread(const B &blocks, const C &centers, const S &screen,
		    const D &density, lockable<F> &fock, Q &queue)
	    : thread_(fock_thread::run<B, C, S, D, F, Q>,
		      blocks, boost::ref(centers), boost::ref(screen),
		      boost::ref(density), boost::ref(fock), boost::ref(queue))
	{
	} 
	void join() { thread_.join(); }

	template<class B, class C, class S, class D, class F, class Q>
	static void run(const B blocks, const C &centers, const S &screen,
			const D &density, lockable<F> &fock, Q &queue) {

	    timer t;

	    std::vector<boost::array<int,4> > quartets;
	    typename Q::value_type task;

	    while (true) {
	    	try { task = queue.backlog().pop(); }
	    	catch(...) {
	    	    try { task = queue.pop(); }
	    	    catch(...) { break; }
	    	}
	    	run_task(blocks, centers, screen, density, fock, task, quartets);
	    }
	    while(true) {
	    	try { task = queue.backlog().pop(); }
	    	catch(...) { break; }
	    	run_task(blocks, centers, screen, density, fock, task, quartets);
	    }

	    std::cout << "CPU thread time: " << t << std::endl;
	}
	template<class B, class C, class S, class D, class F, class T, class Q>
	static void run_task(const B blocks, const C &centers, const S &screen,
			     const D &density, lockable<F> &fock,
			     const T &task, Q &quartets) {
	
	    BOOST_AUTO(ib, &blocks[task[0]]);
	    BOOST_AUTO(kb, &blocks[task[1]]);
	    BOOST_AUTO(jb, &blocks[task[2]]);
	    BOOST_AUTO(lb, &blocks[task[3]]);

	    const BOOST_TYPEOF(ib) block[] = { ib, jb, kb, lb };

	    quartets.clear();
	    screen(block, quartets);
	    if (quartets.size() == 0) return;

	    typedef adapter::rysq::Shell Shell;
	    typedef rysq::Quartet<rysq::Shell> Quartet;
	    Shell a(block[0]->shell());
	    Shell b(block[1]->shell());
	    Shell c(block[2]->shell());
	    Shell d(block[3]->shell());
	    Quartet quartet(a,b,c,d);

	    // if (quartet.L()/2 + 1 < 3) return;
	    // std::cout << quartet << std::endl;

	    int base[] = { ib->firstShell(), jb->firstShell(),
			   kb->firstShell(), lb->firstShell() };

	    rysq::hf::matrix_set<adapter::rysq::density_matrix> D_(base);
	    rysq::hf::matrix_set<adapter::rysq::fock_matrix> F_(base);

	    struct {
		typename F::matrix matrix;
		int index[2];
	    } fock_tls[6];

	    int k = 0;
	    foreach (const size_t (&index)[2],  rysq::index_list()) {
	        size_t i = index[0], j = index[1];
	        size_t im = std::distance(&blocks[0], block[i]);
	        size_t jm = std::distance(&blocks[0], block[j]);

	        D_.get(i,j) = adapter::rysq::density_matrix(density.m(im,jm));
	        D_.get(j,i) = adapter::rysq::density_matrix(density.m(jm,im));

		{
		    typename F::matrix &f = fock.data().m(im,jm);
		    fock_tls[k].matrix.resize(f.size1(), f.size2(),
					      f.block1(), f.block2());
		    fock_tls[k].matrix.clear();
		    fock_tls[k].index[0] = im;
		    fock_tls[k].index[1] = jm;
		    F_.get(i,j) = adapter::rysq::fock_matrix(fock_tls[k].matrix);
		    ++k;
		}
	    }

	    rysq::Parameters parameters(0, -1.0, 4.0, 1e-12);
	    (rysq::Fock(quartet))(centers, quartets, D_, F_, parameters);

	    fock.lock();
	    for (int k = 0; k < 6; ++k) {
		int i = fock_tls[k].index[0], j = fock_tls[k].index[1];
		fock.data().m(i,j) += fock_tls[k].matrix;
	    }
	    fock.unlock();

	}
    };


    struct fock_thread_group {
	static size_t max_threads() {
	    return std::max<size_t>(boost::thread::hardware_concurrency(), 1);
	}
	boost::thread_group group_;
	template<class B, class C, class S, class D, class F, class Q>
	fock_thread_group(const B &blocks, const C &centers, const S &screen,
			  const D &density, lockable<F> &fock, Q &queue,
			  size_t threads = max_threads()) {
	    for (size_t i = 0; i < threads; ++i) {
		group_.add_thread(fock_thread::create(blocks, centers, screen,
						      density, fock, queue));
	    }
	}
	void join_all() { group_.join_all(); }
    };


}

#endif /* _HF_THREAD_HPP_ */
