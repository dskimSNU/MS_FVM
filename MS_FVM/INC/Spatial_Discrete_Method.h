#pragma once
#include <type_traits>
#include <string>

class SDM {};// spatial discrete meethod

namespace ms {
	template <typename T>
	inline constexpr bool is_spatial_discrete_method = std::is_base_of_v<SDM, T>;
}

class FVM : public SDM 
{
public:
	static std::string name(void) { return "FVM"; };
};


class HOM : public SDM {};

// Spatial Discrete Method(SDM)에 변화에 따라 나타나는 Semi_Discrete_Equation(SDE)와 Grid_Info_Extractor(GIE)의 가변성이
// Method 구현 방식의 변화만으로 표현되지 않기 때문에 상속과 합성으로 다형성을 다루지 않고
// SDE와 GIE를 template class로 구현한 뒤 특수화를 통해 다형성을 다룬다.