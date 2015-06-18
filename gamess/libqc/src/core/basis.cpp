#include <functional>
#include <iostream>

#include <boost/typeof/typeof.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
namespace lambda = boost::lambda;

#include "basis.hpp"
#include "util/util.hpp"
#include "iterator/iterator.hpp"


#define OFFSET (-1000)

static std::vector<Basis> bases;

void Basis::add(const Basis::Shell &shell, int center) {
    
    Basis::Shell *shell_ = NULL;

    // search for equivalent shell
    foreach (Basis::Shell *s, this->shells_) {
        if (*s == shell) { shell_ = s; break; }
    }

    // append new shell
    if (!shell_) {
	shells_.push_back(new Shell(shell));
	shell_ = shells_.back();
    }

    size_t size = shell_->size();
    size_t index = this->numShells();
    size_t function = this->numFunctions();

    this->mapped_shells_.push_back(MappedShell(shell_, index, function));
    this->shell_permutation_.push_back(index);
    this->blocks_.push_back(Basis::Block(shell_, index));

    // update function indexes
    std::fill_n(std::back_inserter(this->function_to_shell_), size, index);
    std::generate_n(std::back_inserter(this->function_permutation_),
		    size, util::generator<int>(function));

    this->num_functions_ += size;
    this->max_shell_ = std::max(this->max_shell_, size);

    centers_.push_back((*this->_molecule)(center));

}

template<class C>
typename C::difference_type distance(const C &range,
				     typename C::const_reference r0,
				     typename C::const_reference r1) {
    return (std::find(range.begin(), range.end(), r0) -
	    std::find(range.begin(), range.end(), r1));
}

void Basis::sort() {

    using namespace boost::lambda;

    // sorts shells by size
    std:: sort(shells_.begin(), shells_.end(),
	       bind(&Shell::size, _1) > bind(&Shell::size, _2));

    // sort mapped shells by shell order
    std::sort(mapped_shells_.begin(), mapped_shells_.end(),
	      bind(&distance<BOOST_TYPEOF(shells_)>,
		   boost::ref(shells_), _1, _2) > 0);

    // foreach(MappedShell &mapped, mapped_shells_)
    //     std::cout <<  *mapped.shell << std::endl;

    std::vector<Center> original_centers = this->centers_;

    const Shell* shell = NULL;
    this->blocks_.clear();

    int index = 0;
    int function = 0;
    foreach (MappedShell &mapped, mapped_shells_) {
	mapped.function = function;
	int original_index = mapped.original_index;
	int original_function = mapped.original_function;
	size_t size = mapped.shell->size();
	
	this->shell_permutation_[index] = original_index;
	this->centers_[index] = original_centers[original_index];

        if (mapped.shell != shell) {
            shell = mapped.shell;
            blocks_.push_back(Basis::Block(shell, index, 0));
        }
	++(this->blocks_.back().size_);

        std::fill_n(function_to_shell_.begin() + function, size, index);
        std::generate_n(function_permutation_.begin() + function, size,
			util::generator<int>(original_function));

	++index;
	function += size;
    }

}

std::ostream& operator<<(std::ostream &output, const Basis &b) {
    for (unsigned int i = 0; i < b.shells_.size(); ++i) {
        //output << *(b.shells_.at(i)) << std::endl;
    }

    for (int i = 0; i < b.numShells(); ++i) {
        Basis::MappedShell ms = b.mapped_shells_[i];

        output << i << " -> " << ms.shell;

        output.width(13);
        output.precision(8);
        output << ms.function << ':' << (ms.function + ms.shell->size() - 1) << " ";
        output << ms.original_index << "/" << ms.original_function;

        output.width(13);
        output.precision(8);
        output << b.r(i,0) << "\t" << b.r(i,1) << "\t" << b.r(i,2) << std::endl;
    }

    return output;
}

int Basis_new(int m) {
    Basis b(*Molecule_find_object(m));
    bases.push_back(b);
    return (int(bases.size()) - 1 + OFFSET);
}

Basis* Basis_find_object(int basis) {
    return &bases.at(basis - OFFSET);
}

void Basis_delete(int basis) {
    //bases.erase(bases.begin() + basis + OFFSET);
}

