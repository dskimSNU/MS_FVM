#pragma once
#include "Matrix.h"
#include <string>

class Least_Square
{
public:
	static std::vector<Dynamic_Matrix_> solution_gradients(const std::vector<Dynamic_Matrix_>& solution_delta_matrixes, const std::vector<Dynamic_Matrix_>& least_square_matrixes) {
		const auto num_cell = least_square_matrixes.size();

		std::vector<Dynamic_Matrix_> solution_gradients;
		solution_gradients.reserve(num_cell);
		for (size_t i = 0; i < num_cell; ++i)
			solution_gradients.push_back(solution_delta_matrixes[i] * least_square_matrixes[i]);

		return solution_gradients;
	}
	
	static std::string name(void) { return "Least_Square"; };
};


//template <size_t num_equation, size_t space_domain>
//static std::vector<Matrix<num_equation, space_domain>> solution_gradients(const std::vector<Dynamic_Matrix_>& center_to_center_matrixs, const std::vector<Dynamic_Matrix_>& solution_delta_matrix);

//template <size_t num_equation, size_t space_domain>
//std::vector<Matrix<num_equation, space_domain>> Least_Square::solution_gradients(const std::vector<Dynamic_Matrix_>& center_to_center_matrixes, const std::vector<Dynamic_Matrix_>& solution_delta_matrixes) {
//	const auto num_cell = center_to_center_matrixes.size();
//
//	std::vector<Matrix<num_equation, space_domain>> solution_gradients;
//	solution_gradients.reserve(num_cell);
//	for (size_t i = 0; i < num_cell; ++i) {
//		const auto& Rc = center_to_center_matrixes[i];
//		const auto& dq = solution_delta_matrixes[i];
//		const auto RcT = Rc.transpose();
//		solution_gradients.push_back(dq * RcT * (Rc * RcT).be_inverse());
//	}
//
//	return solution_gradients;
//}
