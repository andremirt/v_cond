#include <map>
#include <utility>
#include <vector>
#include <algorithm>
#include "cxx/sugar.hpp"

namespace rysq {
    namespace memory {

	template<typename T>
	class Allocator {
	public:

	    const size_t page_size;

	    Allocator(size_t page_size) :
		page_size(std::max(page_size/sizeof(T), size_t(1))) { }

	    ~Allocator() {
		foreach (Block *block, blocks) delete block;
	    }

	    T* push(size_t n) {
		if (blocks.empty()) blocks.push_back(new Block(n, page_size));
	
		Block *block = blocks.back();
		if (block->current + n >= block->last) {
		    blocks.push_back(new Block(n, page_size));
		    block = blocks.back();
		}

		T *ptr = block->current;
		block->current +=n;
		return ptr;
	    }


	private:
	    class Block {
	    public:
		T *ptr, *last, *current;
		Block(size_t n, size_t page_size) {
		    size_t size = std::max(n, page_size);
		    current = ptr = new T[size];
		    last = ptr + size;
		    std::fill( ptr, last, 0);
		}
		~Block() { delete ptr; }

	    };


	    std::vector<Block*> blocks;

	};


	template <typename T>
	class Map {

	public:

	    typedef std::pair<T*,size_t> Value;
	    typedef std::pair<uint64_t, Value> KeyValue;
	    typedef std::map<uint64_t, Value> map;
	    typedef typename map::iterator Iterator;

	    struct Entry {
		int i,j;
		T *ptr;
		size_t n;

		Entry(KeyValue &kv)
		    : i(int(kv.first & 0xffffffff)),
		      j(int(kv.first >> 32 & 0xffffffff)),
		      ptr(kv.second.first),
		    n(kv.second.second) {
		}
	    };

	    Map( size_t page_size = 8192) : allocator(page_size) {}

	    uint64_t key(uint32_t i, uint32_t j) {
		return uint64_t(i) | uint64_t(j) << 32;
	    }

	    T* get(uint32_t i, uint32_t j) {
		Iterator it = map_.find(key(i,j));
		if (it != map_.end()) return ((*it).second).first;
		return NULL;
	    }

	    T* get(uint32_t i, uint32_t j, size_t n) {
		T *ptr = get(i,j);
		if (!ptr) {
		    assert(map_.size() < map_.max_size());
		    ptr = allocator.push(n);
		    map_[key(i,j)] = Value(ptr,n);
		}
		return ptr;
	    }


	    Entry pop() {
		Iterator it = map_.begin();
		KeyValue kv = *it;
		map_.erase(kv.first);
		return kv;
	    }

	    bool empty() { return map_.empty(); }

	private:

	    map map_;
	    Allocator<T> allocator;

	};

    }
}
