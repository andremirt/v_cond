#ifndef _BASIS_H
#define _BASIS_H

/**
   @file
   @brief
*/

#include "molecule.hpp"
#include "util/numeric.hpp"

#include <string>
#include <vector>
#include <iostream>
#include <math.h>

#include "core/shell.hpp"
#include "core/block.hpp"

/**
   @brief
*/
class Basis {

    /**
       @brief
    */
    friend std::ostream& operator<<(std::ostream &output, const Basis &b);

public:

    typedef basis::Shell Shell;
    typedef basis::Block Block;
    typedef Molecule::Center Center;

    /**
       @brief Basis constructor
       @param molecule
    */
    Basis(const Molecule &molecule):
	_molecule(&molecule), max_shell_(0), num_functions_(0), _sorted(false) {}

    /** @brief Get number of shells */
    int numShells() const { return mapped_shells_.size(); }

    /** @brief Get total number of functions */
    int numFunctions() const { return num_functions_; }

    /** @brief Get the shells */
    const Shell& operator()(int shell) const { return *mapped_shells_[shell].shell; }

    /**
       @brief
       @param K
       @param exp
       @param C
       @param L
       @param center
    */
    void add(int K, double *exp, double *C, int L, int center) {
	add(Shell(K, exp, C, L), center);
    }

    /**
       @brief
       @param K
       @param exp
       @param Cs
       @param Cp
       @param center
    */
    void add(int K, double *exp, double *Cs, double *Cp, int center) {
    add(Shell(K, exp, Cs, Cp), center);
    }

    /**
       @brief
       @param i
    */
    int firstShellFunction(int i) const { return mapped_shells_[i].function; }
    /**
       @brief
     */
    int maxShell() const { return max_shell_; }
    /**
       @brief
       @param shell
       @param coord
    */
    double r(int shell, int coord) const { return this->r(shell)[coord]; }
    /**
       @brief
       @param shell
     */
    const double* r(int shell) const { return &centers_.at(shell)[0]; }

    const std::vector<Center>& centers() const { return centers_; }

    /**
       @brief
       @param shell
    */
    //int functionToShell(int shell) const { return function_to_shell_[shell]; }

    void sort();

    int shell_permutation(int i) const { return shell_permutation_.at(i); }

    const std::vector<int>& shell_permutations() const {
	return shell_permutation_;
    }

    const std::vector<int>& function_permutations() const {
	return function_permutation_;
    }

    int function_permutation(int i) const { return function_permutation_.at(i); }

    std::vector<int> getShellDimensions() const {
        std::vector<int> blockDims;

        for (int i = 0; i < numShells(); ++i) {
            blockDims.push_back((*this)(i).size());
        }

        return blockDims;
    }

     int maxBlock() const {
	return (std::max_element(blocks_.begin(), blocks_.end(),
				 Basis::Block::compareSize))->numShells();
    }

//     const std::vector<int>& shellPermutation() const {
// 	return shellPermutation_;
//     }

    typedef std::vector<Block>::iterator block_iterator;
    typedef std::vector<Block>::const_iterator const_block_iterator;

    const Block& block(int i) const { return blocks_.at(i); }
    const std::vector<Block>& blocks() const { return blocks_; }
    const Block* firstBlock() const { return &(blocks_.front()); }
    const Block* lastBlock() const { return &(blocks_.back()); }

    bool isSorted() const { return _sorted; }

    template<class T> void normalize(T &A) const;

    /** order shell such that all shells of the same type back contiguous 
	does this to _shells and _shellCenters
	@brief
	@param output
	@param m
	@param a
	@param b
	@param shell
	@param center
	@param firstFunction
	@param firstFunction0
    */

private:
    struct MappedShell {
	Shell *shell; //Change to index to vector
	int original_index;
	int original_function; //Const
	int function;

	/**
	   @brief
	   @param shell
	   @param center
	   @param index
	*/
        MappedShell(Shell *shell, int index, int function)
	    : shell(shell),
	      original_index(index),
	      original_function(function),
	      function(function) {}
	
	typedef Shell* pointer;
	operator pointer() const { return shell; }

	bool operator>(const MappedShell &rhs) const {
	    return (this->shell > rhs.shell);
	}

	bool operator>=(const MappedShell &rhs) const {
	    return (this->shell->L() >= rhs.shell->L());
	}

    };


    /** creates blocks based contiguous _shells of same type in same block
	@brief
	@param s
	@param center
    */
    //void blockShells();

    void add(const Shell &s, int center);

    /** Molecule Object
     @brief
     @param molecule
     */
    const Molecule *_molecule;
    /** @brief List of shells*/
    std::vector<Shell*> shells_;

    /** @brief List of mapped shells*/
    std::vector<MappedShell> mapped_shells_;

    std::vector<Center> centers_;

    std::vector<int> function_to_shell_;
    std::vector<int> function_permutation_;
    std::vector<int> shell_permutation_;

    std::vector<Block> blocks_;

    size_t max_shell_;

    size_t num_functions_;

    bool _sorted;
};




extern "C" {

    /**
       @brief
       @param m
    */
    int Basis_new(int m);

    /**
       @brief
       @param basis
    */
    Basis* Basis_find_object(int basis);

    /**
       @brief
       @param
     */
    void Basis_delete(int basis);
}



#endif // _BASIS_H
