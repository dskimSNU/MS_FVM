#include "../INC/Governing_Equation.h"

Linear_Advection_2D::Physical_Flux_ Linear_Advection_2D::physical_flux(const Solution_& solution) {
	const auto [x_advection_speed, y_advection_speed] = Linear_Advection_2D::advection_speeds_;
	const auto sol = solution[0];

	return { x_advection_speed * sol , y_advection_speed * sol };
}

std::vector<Linear_Advection_2D::Physical_Flux_> Linear_Advection_2D::physical_fluxes(const std::vector<Solution_>& solutions) {
	//static size_t num_solution = solutions.size();
	const size_t num_solution = solutions.size();

	const auto [x_advection_speed, y_advection_speed] = Linear_Advection_2D::advection_speeds_;
	std::vector<Physical_Flux_> physical_fluxes(num_solution);
	for (size_t i = 0; i < num_solution; ++i) {
		const auto sol = solutions[i][0];
		physical_fluxes[i] = { x_advection_speed * sol , y_advection_speed * sol };
	}

	return physical_fluxes;
}

std::vector<std::array<double, Linear_Advection_2D::space_dimension_>> Linear_Advection_2D::coordinate_projected_maximum_lambdas(const std::vector<Solution_>& solutions) {
	//static size_t num_solution = solutions.size();
	const size_t num_solution = solutions.size();

	static double absolute_x_advection_speed = std::abs(advection_speeds_[0]);
	static double absolute_y_advection_speed = std::abs(advection_speeds_[1]);

	std::vector<std::array<double, Linear_Advection_2D::space_dimension_>> projected_maximum_lambdas(num_solution, { absolute_x_advection_speed,absolute_y_advection_speed });
	return projected_maximum_lambdas;
}

double Linear_Advection_2D::inner_face_maximum_lambda(const Solution_& solution_o, const Solution_& solution_n, const Space_Vector_& nomal_vector) {
	return std::abs(nomal_vector.inner_product(advection_speeds_));
}


Burgers_2D::Physical_Flux_ Burgers_2D::physical_flux(const Solution_& solution) {
	const auto sol = solution[0];

	const auto temp_val = 0.5 * sol * sol;
	return { temp_val * temp_val };
}

std::vector<Burgers_2D::Physical_Flux_> Burgers_2D::physical_fluxes(const std::vector<Solution_>& solutions) {
	static size_t num_solution = solutions.size();


	std::vector<Physical_Flux_> physical_fluxes(num_solution);
	for (size_t i = 0; i < num_solution; ++i) {
		const auto sol = solutions[i][0];
		const auto temp_val = 0.5 * sol * sol;
		physical_fluxes[i] = { temp_val, temp_val };
	}

	return physical_fluxes;
}

std::vector<std::array<double, Burgers_2D::space_dimension_>> Burgers_2D::coordinate_projected_maximum_lambdas(const std::vector<Solution_>& solutions) {
	static size_t num_solution = solutions.size();

	std::vector<std::array<double, Burgers_2D::space_dimension_>> projected_maximum_lambdas(num_solution);
	for (size_t i = 0; i < num_solution; ++i) {
		const auto maximum_lambdas = std::abs(solutions[i][0]);
		projected_maximum_lambdas[i] = { maximum_lambdas, maximum_lambdas };
	}

	return projected_maximum_lambdas;
}

double Burgers_2D::inner_face_maximum_lambda(const Solution_& solution_o, const Solution_& solution_n, const Space_Vector_& nomal_vector) {
	const auto normal_component_sum = nomal_vector[0] + nomal_vector[1];	
	return std::max(std::abs(solution_o[0] * normal_component_sum), std::abs(solution_n[0] * normal_component_sum));
}