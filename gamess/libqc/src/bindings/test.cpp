#include <boost/bind.hpp>
#include <algorithm>
#include <string>
#include <vector>

int test() {
    std::vector<std::string> vec;
    std::less<size_t> cmp;
    std::remove_if(vec.begin(),
		   vec.end(),
		   boost::bind(cmp, 5, boost::bind(&std::string::length, _1)));
}
