#pragma once
#include <type_traits>

template <double value>
class TSM { //time step method
public:
	static constexpr double constant() { return value; };
}; 

namespace ms {
	template <typename T>
	inline constexpr bool is_time_step_method = std::is_base_of_v<TSM<T::constant()>, T>;
}


template <double value>
class CFL : public TSM<value> {};


template <double value>
class ConstantDt : public TSM<value> {};