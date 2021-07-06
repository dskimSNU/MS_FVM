#include "../INC/Governing_Equation.h"

std::vector<Scalar_Conservation_Law::Physical_Flux> Linear_Advection_2D::calculate_physical_fluxes(const std::vector<Solution>& solutions) {
	static size_t num_solution = solutions.size();
	
	std::vector<Physical_Flux> physical_fluxes(num_solution);
	for (size_t i = 0; i < num_solution; ++i) {
		const auto solution = solutions[i][0];
		physical_fluxes[i] = { advection_speeds_[0] * solution, advection_speeds_[1] * solution };
	}

	return physical_fluxes;
}

std::array<std::vector<double>, 2> Linear_Advection_2D::calculate_coordinate_projected_maximum_lambdas(const std::vector<Solution>& solutions) {
	static size_t num_solution = solutions.size();
	static double absolute_x_advection_speed = std::abs(advection_speeds_[0]);
	static double absolute_y_advection_speed = std::abs(advection_speeds_[1]);

	std::vector<double> x_projected_maximum_lambdas(num_solution, absolute_x_advection_speed);
	std::vector<double> y_projected_maximum_lambdas(num_solution, absolute_y_advection_speed);
	return { x_projected_maximum_lambdas, y_projected_maximum_lambdas };
}

double Linear_Advection_2D::calculate_inner_face_maximum_lambdas(const Solution& solution_o, const Solution& solution_n, const Physical_Domain_Vector& nomal_vector) {
	return std::abs(nomal_vector.inner_product(advection_speeds_));
}