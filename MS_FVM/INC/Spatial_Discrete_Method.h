#pragma once
#include <type_traits>

//dummy class for specialize cells and boundaries and innerfaces RHS calculation algorithm.
class Spatial_Discrete_Method {};

namespace ms {
	template <typename T>
	inline constexpr bool is_spatial_discrete_method = std::is_base_of_v<Spatial_Discrete_Method, T>;
}

class FVM : public Spatial_Discrete_Method {};
class HOM : public Spatial_Discrete_Method {};

