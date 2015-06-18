#ifndef RUNTIME_HPP
#define RUNTIME_HPP

#include <string>
#include <vector>
#include <cstdlib>
#include <boost/lexical_cast.hpp>

namespace runtime {

    template<typename T>
    T num_threads(const T &default_value, const char *variable = "NUM_THREADS") {
	char *value = std::getenv(variable);
	if (value) return boost::lexical_cast<T>(value);
	else return default_value;
    }

    template<typename T>
    std::vector<T> cuda_devices(const std::vector<T> &default_value,
				const char *variable = "CUDA_DEVICES") {
	if (!std::getenv(variable)) return default_value;

	std::string value(std::getenv(variable));
	std::vector<T> v;
	while (!value.empty()) {
	    size_t token = value.find(",");
	    T t = boost::lexical_cast<T>(value.substr(0,token));
	    if (t < T(0)) return std::vector<T>();
	    v.push_back(t);
	    token -= (token == std::string::npos);
	    value.erase(0, token+1);
	}
	return v;
    }

}

#endif /* RUNTIME_HPP */
