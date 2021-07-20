#pragma once
#include <type_traits>
#include <string>

#include "Gradient_Method.h"

class RM {};	// Reconstruction Method

class Constant_Reconstruction : public RM 
{
public:
	static std::string name(void) { return "Constant_Reconstruction"; };
};


template <typename Gradient_Method>
class Linear_Reconstruction : public RM
{
public:
	static std::string name(void) { return "Linear_Reconstruction_" + Gradient_Method::name(); };
};


class MLP : public RM {};

template <typename Gradient_Method>
class MLP_u1 : public MLP 
{
public:
	static std::string name(void) { return "MLP_u1_" + Gradient_Method::name(); };
};

namespace ms {
	template <typename T>
	inline constexpr bool is_reconsturction_method = std::is_base_of_v<RM, T>;

	template <typename T>
	inline constexpr bool is_MLP_method = std::is_base_of_v<MLP, T>;
}


