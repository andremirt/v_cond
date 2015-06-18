#include <map>
#include <vector>
#include <algorithm>
#include <iterator>
#include <boost/array.hpp>

#include "cxx/sugar/sugar.hpp"
#include "cxx/iterator/product.hpp"

#include "cudaxx/cudaxx.hpp"

#include <cuda.h>
#include <cuda_runtime.h>

#include <rysq.hpp>

using namespace cudaxx;

// const device::const_array<ushort3> index(const rysq::Quartet< rysq:: Shell> &quartet) {
// }

class Index {

public:
    static Index* instance() {
	if (instance_ == NULL)  {
	    std::vector<rysq::type> types;
	    for (int i = 0; i <= RYSQ_LMAX; ++i) {
		types.push_back(rysq::type(i));
	    }
	    types.push_back(rysq::SP);
	    instance_ = new Index(types);
	}
	return instance_;
    }

private:
    typedef boost::array<rysq::type,4> type_quartet;
	typedef std::vector<ushort3> index_vector;
    std::map<type_quartet, size_t> index_map;
    device_ptr<ushort3> index_;
    static Index* instance_;

    Index(const std::vector<rysq::type> &types) {

	index_vector index(1024);
	index_vector::iterator it = index.begin();

	std::ostream_iterator<ushort3, char> ushort3_output(std::cout, "\n");
	std::ostream_iterator< int > key_output(std::cout, " " );

	foreach (const  type_quartet q,  cxx::iterator::product<4>(types)) {
	    size_t size =  rysq::Quartet<rysq::Shell>::size(q.elems);
	    if (size > 2048) continue;

	    if (index_map.find(q) == index_map.end()) {
		index_vector qindex;
	    	// initialize(key, index);
	    	assert(index_map.size() < index_map.max_size());
	    	index_map[q] = index.size();
		std::copy(qindex.begin(), qindex.end(), std::back_inserter(index));
	    }
	}
	 // index_ = indexv;

    }

    // const ushort3* get(const Key &key) const {
    // 	std::map<Key, size_t>::const_iterator it = index_map.find(key);
    // 	if (it == index_map.end()) return NULL;
    // 	return gIndexes->data().begin() + (*it).second;
    // }

    // static int indexOf(boost::array<int,4> idx, boost::array<int,4> L) {
    // 	int m0 = L[0]+1;
    // 	int m1 = L[1]+1;
    // 	int m2 = L[2]+1;
    // 	return idx[0] + (idx[1] + (idx[2] + idx[3]*m2)*m1)*m0;
    // }

    static index_vector initialize(const type_quartet &quartet) {
    	index_vector index(256);
    	// foreach (const type_quartet q, product(quartet)) {
    	//     int xyz4[3][4];
    	//     for (XRANGE(i, 4)) {
    	// 	xyz4[0][i] = rysq::LX[q[i]];
    	// 	xyz4[1][i] = rysq::LY[q[i]];
    	// 	xyz4[2][i] = rysq::LZ[q[i]];
    	//     }
    	//     index.push_back(index3(xyz4, q));
    	// }
    	return index;
    }

};
