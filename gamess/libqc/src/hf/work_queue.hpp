#ifndef _HF_WORK_QUEUE_HPP_
#define _HF_WORK_QUEUE_HPP_

#include <queue>
#include <boost/array.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/locks.hpp>
#include <boost/thread/condition_variable.hpp>
#include <boost/mpl/vector_c.hpp>

#include "generator.hpp"
#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <algorithm>
#include "externals/cxx/array.hpp"

namespace hf {

    struct work_queue {
	typedef boost::array<int,4> value_type;
	int N;
	value_type current;
	typedef generator::For<int,4> generator_type;
	generator_type g;
	
	static generator_type make_generator(size_t N) {
	    using generator::range;
	    using generator::make_for;
	    using namespace boost::lambda;
	    const int&(*max)(const int&, const int&) = std::max;
	    return generator_type(N, N, range(bind(max, _1, _2), N), range(_2, _3+1),
				  boost::mpl::vector_c<size_t,3,2,1,0>());
	}
	work_queue(size_t N,size_t throttle = 8196)
	    : g(make_generator(N)), backlog_(throttle)
	{
	    this->N = N;
	    current.assign(0);
	}
	bool empty() const { return current[3] >= N; }
	value_type front() {
	    //boost::mutex::scoped_lock(mutex_);
	    boost::lock_guard<boost::mutex> lock(mutex_);
	    return current;
	}
	value_type pop() {

	    //boost::mutex::scoped_lock(mutex_);
	    boost::lock_guard<boost::mutex> lock(mutex_);
	    // std::cout << g.next() << current << std::endl;
	    if (empty()) throw std::exception();

	    value_type next = current;

	    // for (RANGE(lb, basis.blocks())) {
	    // 	for (RANGE(jb, basis.blocks())) {
	    // 	    for (RANGE(kb, std::max(lb, jb), basis. blocks().end())) {
	    // 		for (RANGE(ib, jb, kb+1)) {

	    bool advance;
	    current[0] += 1;

	    advance = (current[0] >= (current[1]+1));
	    if (advance) {
		current[1] += 1;
		current[0] = current[2];
	    }

	    advance = advance && (current[1] == N);
	    if (advance) {
		current[2] += 1;
		current[1] = std::max(current[2], current[3]);
		current[0] = current[2];
	    }

	    advance = advance && (current[2] == N);
	    if (advance) {
		current[3] += 1;
		current[2] = 0;
		current[1] = std::max(current[2], current[3]);
		current[0] = current[2];
	    }

	    //std::cout << next << std::endl;
	    return next;
	}
	struct backlog_queue {
	    explicit backlog_queue(size_t throttle = 8196) {
		throttle_ = std::max<size_t>(throttle, 2);
		throttle_set_ = false;
	    }
	    bool empty () const { return queue_.empty(); }
	    void push(const value_type &work) {
		// boost::mutex::scoped_lock throttle_lock(throttle_mutex_);
		boost::mutex::scoped_lock lock(mutex_);
		while (throttle_set_) throttle_condition_.wait(lock);
		queue_.push(work);
		if (queue_.size() > throttle_) {
		    //std::cout << "throttle set" << std::endl;
		    throttle_set_ = true;
		}
	    }
	    value_type pop() {
		//boost::unique_lock<boost::mutex> lock(mutex_);
		 boost::mutex::scoped_lock lock(mutex_);
		if (queue_.empty()) {
		    throw std::exception();
		}
		value_type front = queue_.front();
		queue_.pop();
		if (queue_.size() < throttle_/2) {
		    throttle_set_ = false;
		    throttle_condition_.notify_all();
		}
		// std::cout <<  queue_.size()<< front << std::endl;
		return front;
	    }
	private:
	    std:: queue< value_type> queue_;
	    size_t throttle_;
	    boost::mutex mutex_;
	    bool throttle_set_;
	    boost::condition_variable throttle_condition_;
	};
	backlog_queue& backlog() { return backlog_; }
    private:
	boost::mutex mutex_;
	backlog_queue backlog_;
    };

}


// for (RANGE(lb, N)) {
// 	for (RANGE(jb, N)) {
// 	    for (RANGE(kb, std::max(lb, jb), N)) {
// 		for (RANGE(ib, jb, kb+1)) {



#endif /* _HF_WORK_QUEUE_HPP_ */
