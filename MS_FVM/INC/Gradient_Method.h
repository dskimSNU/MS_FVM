#pragma once
#include "Matrix.h"
#include <string>

class Least_Square
{
public:
	static std::vector<Dynamic_Matrix_> solution_gradients(const std::vector<Dynamic_Matrix_>& center_to_center_matrixs, const std::vector<Dynamic_Matrix_>& solution_delta_matrix);
	static std::string name(void) { return "Least_Square"; };
};