#pragma once
#include <type_traits>

#include "Gradient_Method.h"

class RM {};	// Reconstruction Method

class Constant_Reconstruction : public RM {};

class MLP : public RM {};

template <typename Gradient_Method>
class MLP_u1 : public MLP {};

namespace ms {
	template <typename T>
	inline constexpr bool is_reconsturction_method = std::is_base_of_v<RM, T>;

	template <typename T>
	inline constexpr bool is_MLP_method = std::is_base_of_v<MLP, T>;
}


