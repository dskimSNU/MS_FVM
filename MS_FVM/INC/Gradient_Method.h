#pragma once
#include "Matrix.h"

class Least_Square
{
public:
	static std::vector<Dynamic_Matrix_> solution_gradients(const std::vector<Dynamic_Matrix_>& center_to_center_matrixs, const std::vector<Dynamic_Matrix_>& solution_delta_matrix);
};