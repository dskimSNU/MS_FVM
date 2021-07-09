#include "../INC/Governing_Equation.h"

std::vector<Scalar_Conservation_Law_2D::Physical_Flux> Linear_Advection_2D::calculate_physical_fluxes(const std::vector<Solution>& solutions) {
	static size_t num_solution = solutions.size();

	std::vector<Physical_Flux> physical_fluxes(num_solution);
	for (size_t i = 0; i < num_solution; ++i) {
		const auto solution = solutions[i][0];
		physical_fluxes[i] = { advection_speeds_[0] * solution, advection_speeds_[1] * solution };
	}

	return physical_fluxes;
}

std::vector<std::array<double, Linear_Advection_2D::physical_domain_dimension_>> Linear_Advection_2D::calculate_coordinate_projected_maximum_lambdas(const std::vector<Solution>& solutions) {
	static size_t num_solution = solutions.size();

	static double absolute_x_advection_speed = std::abs(advection_speeds_[0]);
	static double absolute_y_advection_speed = std::abs(advection_speeds_[1]);

	std::vector<std::array<double, Linear_Advection_2D::physical_domain_dimension_>> projected_maximum_lambdas(num_solution, { absolute_x_advection_speed,absolute_y_advection_speed });
	return projected_maximum_lambdas;
}

double Linear_Advection_2D::calculate_inner_face_maximum_lambda(const Solution& solution_o, const Solution& solution_n, const Physical_Domain_Vector& nomal_vector) {
	return std::abs(nomal_vector.inner_product(advection_speeds_));
}

std::vector<Scalar_Conservation_Law_2D::Physical_Flux> Burgers_2D::calculate_physical_fluxes(const std::vector<Solution>& solutions) {
	static size_t num_solution = solutions.size();

	std::vector<Physical_Flux> physical_fluxes(num_solution);
	for (size_t i = 0; i < num_solution; ++i) {
		const auto temp_val = 0.5 * solutions[i][0] * solutions[i][0];
		physical_fluxes[i] = { temp_val, temp_val };
	}

	return physical_fluxes;
}

std::vector<std::array<double, Burgers_2D::physical_domain_dimension_>> Burgers_2D::calculate_coordinate_projected_maximum_lambdas(const std::vector<Solution>& solutions) {
	static size_t num_solution = solutions.size();

	std::vector<std::array<double, Burgers_2D::physical_domain_dimension_>> projected_maximum_lambdas(num_solution);
	for (size_t i = 0; i < num_solution; ++i) {
		const auto maximum_lambdas = std::abs(solutions[i][0]);
		projected_maximum_lambdas[i] = { maximum_lambdas, maximum_lambdas };
	}

	return projected_maximum_lambdas;
}

double Burgers_2D::calculate_inner_face_maximum_lambda(const Solution& solution_o, const Solution& solution_n, const Physical_Domain_Vector& nomal_vector) {
	const auto normal_component_sum = nomal_vector[0] + nomal_vector[1];	
	return std::max(std::abs(solution_o[0] * normal_component_sum), std::abs(solution_n[0] * normal_component_sum));
}