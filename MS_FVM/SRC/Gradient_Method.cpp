#include "../INC/Gradient_Method.h"

std::vector<Dynamic_Matrix_> Least_Square::solution_gradients(const std::vector<Dynamic_Matrix_>& center_to_center_matrixes, const std::vector<Dynamic_Matrix_>& solution_delta_matrixes) {
	const auto num_cell = center_to_center_matrixes.size();

	std::vector<Dynamic_Matrix_> solution_gradients;
	solution_gradients.reserve(num_cell);
	for (size_t i=0; i<num_cell; ++i){
		const auto& Rc = center_to_center_matrixes[i];
		const auto& dQ = solution_delta_matrixes[i];
		const auto RcT = Rc.transpose();
		const auto pseudo_inverse = RcT * (Rc * RcT).be_inverse();


		solution_gradients.push_back(dQ * pseudo_inverse);
	}

	return solution_gradients;
}
