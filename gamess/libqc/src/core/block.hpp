#ifndef _BASIS_BLOCK_HPP_
#define _BASIS_BLOCK_HPP_

#include "core/shell.hpp"

namespace basis {

    class Block {
	friend class Basis;

	const Shell *shell_;
	int index_;


    public:
	int size_;
	Block(const Shell *shell, int firstShell, int numShells = 1)
	    : shell_(shell),
	      index_(firstShell),
	      size_(numShells) {
	}

	size_t size() const { return size_; }

	const Shell& shell() const { return *shell_; }

	int numShells() const { return size_; }


	int firstShell() const { return index_; }

	int lastShell() const { return index_ + size_ - 1; }

	int numFunctions() const { return size_*shell_->size(); }

	static bool compareSize(const Block &a, const Block &b) {
	    return a.numShells() > b.numShells();
	}

	bool operator==(const Block &b) const { return (this == &b); }

    };

}

#endif /* _BASIS_BLOCK_HPP_ */
