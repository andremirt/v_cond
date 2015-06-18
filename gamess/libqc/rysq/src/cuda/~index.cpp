#ifndef _RYSQ_CUDA_KERNELS_INDEX_H_
#define _RYSQ_CUDA_KERNELS_INDEX_H_

#include <map>
#include <vector>
#include <boost/array.hpp>
#include "sugar/sugar.hpp"

namespace rysq {
    namespace cuda {

// 	std::ostream& operator<< (std::ostream& os, ushort3& v) {
// 	    return os << v.x << ", "<< v.y << ", " << v.z;
// 	}

	class Index {

	public:
	    static Index* instance() {
		if (instance_ == NULL)  {
		    std::vector<int> types;
		    for (int i = 0; i <= RYSQ_LMAX; ++i) {
		        types.push_back(i);
		    }
		    types.push_back(-1);
		    instance_ = new Index(types);
		}
		return instance_;
	    }

	    static const ushort3* get(int i, int j, int k = 0, int l = 0) {
		Key key = {{ i, j, k, l }};
		return instance()->get(key);
	    }

	private:
	    typedef boost::array<int,4> Key;
	    std::map<Key, size_t> indexMap;
	    std::vector<ushort3> indexes;
	    device_array<ushort3> *gIndexes;
	    static Index* instance_;

	    Index(const std::vector<int> &types) {
		indexes.reserve(1024);
		std::vector<ushort3>::iterator it = indexes.begin();

		std::ostream_iterator<ushort3, char> ushort3_output(std::cout, "\n");
		std::ostream_iterator< int > key_output(std::cout, " " );

		foreach(int t0, types) {
		    // size_t size0 = shell::size(t0);
		    foreach(int t1, types) {
			// size_t size1 = shell::size(t1);
			foreach(int t2, types) {
			    // size_t size2 = shell::size(t2);
			    foreach(int t3, types) {
				// size_t size3 = shell::size(t3);

				size_t size = Quartet<Shell>::size(t0, t1, t2, t3);
				if (size > 2048) continue;

				Key key = {{ t0, t1, t2, t3 }};
				if (indexMap.find(key) == indexMap.end()) {
				    ushort3 *index  = new ushort3[size];
				    initialize(key, index, size);
				    assert(indexMap.size() < indexMap.max_size());
				    indexMap[key] = indexes.size();
				    indexes.insert(indexes.end(), index, index + size);

// 				    std::cout << "key: ";
// 				    std::copy(key.elems, key.elems + 4, key_output);
// 				    std::cout << "\n";
// 				    for (int i = 0; i < size; ++i) {
// 					std::cout << index[i] << "\n";
// 				    }

				    delete[] index;
				    
				}

			    }
			}
		    }
		}
		gIndexes = new device_array<ushort3>(indexes);

	    }

	    const ushort3* get(const Key &key) const {
		std::map<Key, size_t>::const_iterator it = indexMap.find(key);
		if (it == indexMap.end()) return NULL;
		return gIndexes->data().begin() + (*it).second;
	    }

	    static int indexOf(boost::array<int,4> idx, boost::array<int,4> L) {
		int m0 = L[0]+1;
		int m1 = L[1]+1;
		int m2 = L[2]+1;
		return idx[0] + (idx[1] + (idx[2] + idx[3]*m2)*m1)*m0;
	    }

	    static void initialize(Key type, ushort3 *index, size_t size) {

		boost::array<int,4> L = {{ abs(type[0]), abs(type[1]),
					   abs(type[2]), abs(type[3]) }};

		shell a(type[0]), b(type[1]), c(type[2]), d(type[3]);

		size_t ijkl = 0;
		for (RANGE(l, d)) {
		    int lx = rysq::LX[l];
		    int ly = rysq::LY[l];
		    int lz = rysq::LZ[l];

		    for (RANGE(k, c)) {
			int kx = rysq::LX[k];
			int ky = rysq::LY[k];
			int kz = rysq::LZ[k];

			for (RANGE(j, b)) {
			    int jx = rysq::LX[j];
			    int jy = rysq::LY[j];
			    int jz = rysq::LZ[j];

			    for (RANGE(i, a)) {
				int ix = rysq::LX[i];
				int iy = rysq::LY[i];
				int iz = rysq::LZ[i];

				boost::array<int,4> x = {{ ix, jx, kx, lx }};
				boost::array<int,4> y = {{ iy, jy, ky, ly }};
				boost::array<int,4> z = {{ iz, jz, kz, lz }};

				assert(ijkl < size);

				index[ijkl].x = Index::indexOf(x, L);
				index[ijkl].y = Index::indexOf(y, L);
				index[ijkl].z = Index::indexOf(z, L);
				++ijkl;
			    }
			}
		    }
		}
		assert(ijkl == size);
	    }

	};

	cuda::Index* cuda::Index::instance_ = NULL;

    }
}
 
#endif /* _RYSQ_CUDA_KERNELS_INDEX_H_ */
