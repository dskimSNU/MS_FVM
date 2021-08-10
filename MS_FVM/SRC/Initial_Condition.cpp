#include "../INC/Inital_Condition.h"

std::vector<Sine_Wave_2D::Solution> Sine_Wave_2D::calculate_solutions(const std::vector<Space_Vector_>& cell_centers) {
	const auto num_cell = cell_centers.size();
	
	std::vector<Solution> solutions_(num_cell);
	for (size_t i = 0; i < num_cell; ++i) {
		const auto x_coord = cell_centers[i][0];
		const auto y_coord = cell_centers[i][1];
		//solutions_[i] = std::sin(2 * pi_ * x_coord) * std::sin(2 * pi_ * y_coord);

		const auto x_wave_number = 2 * pi_ / x_wave_length_;
		const auto y_wave_number = 2 * pi_ / y_wave_length_;

		solutions_[i] = std::sin(x_wave_number * x_coord ) * std::sin(y_wave_number * y_coord);
	}

	return solutions_;
}

std::string Sine_Wave_2D::name(void) {
	return "Sine_Wave_2D_WL(" + ms::d_to_str_with_precision(x_wave_length_, 2) + "," + ms::d_to_str_with_precision(y_wave_length_, 2) + ")"; 
};


template <>
std::vector<Sine_Wave_2D::Solution> Sine_Wave_2D::calculate_exact_solutions<Linear_Advection_2D>(const std::vector<Space_Vector_>& cell_centers, const double end_time){
	const auto num_cell = cell_centers.size();
	const auto [x_advection_speed, y_advection_speed] = Linear_Advection_2D::advection_speed();

	std::vector<Solution> exact_solutions_(num_cell);
	for (size_t i = 0; i < num_cell; ++i) {
		const auto x_coord = cell_centers[i][0];
		const auto y_coord = cell_centers[i][1];
		exact_solutions_[i] = std::sin(2 * pi_ * (x_coord - x_advection_speed * end_time)) * std::sin(2 * pi_ * (y_coord - y_advection_speed * end_time));
	}

	return exact_solutions_;
}


std::vector<Square_Wave_2D::Solution> Square_Wave_2D::calculate_solutions(const std::vector<Space_Vector_>& cell_centers) {
	const auto num_cell = cell_centers.size();

	std::vector<Solution> solutions(num_cell);
	for (size_t i = 0; i < num_cell; ++i) {
		const auto x_coord = cell_centers[i][0];
		const auto y_coord = cell_centers[i][1];

		if (0.25 <= x_coord && x_coord <= 0.75 && 0.25 <= y_coord && y_coord <= 0.75)
			solutions[i] = 1;
		else
			solutions[i] = 0;
	}

	return solutions;
}


template <>
std::vector<Square_Wave_2D::Solution> Square_Wave_2D::calculate_exact_solutions<Linear_Advection_2D>(const std::vector<Space_Vector_>& cell_centers, const double end_time) {
	const auto num_cell = cell_centers.size();
	const auto [x_advection_speed, y_advection_speed] = Linear_Advection_2D::advection_speed();

	std::vector<Solution> exact_solutions_(num_cell);
	for (size_t i = 0; i < num_cell; ++i) {
		const auto x_coord = cell_centers[i][0];
		const auto y_coord = cell_centers[i][1];

		//Assume that domian [0,1] x [0,1]
		const auto exact_x_start	= 0.25 + x_advection_speed * end_time - static_cast<int>(0.25 + x_advection_speed * end_time);
		const auto exact_x_end		= 0.75 + x_advection_speed * end_time - static_cast<int>(0.75 + x_advection_speed * end_time);
		const auto exact_y_start	= 0.25 + y_advection_speed * end_time - static_cast<int>(0.25 + y_advection_speed * end_time);
		const auto exact_y_end		= 0.75 + y_advection_speed * end_time - static_cast<int>(0.75 + y_advection_speed * end_time);

		if (exact_x_start <= x_coord && x_coord <= exact_x_end &&
			exact_y_start <= y_coord && y_coord <= exact_y_end)
			exact_solutions_[i] = 1;
		else
			exact_solutions_[i] = 0;
	}

	return exact_solutions_;
}


std::vector<Modifid_SOD_2D::Solution_> Modifid_SOD_2D::calculate_solutions(const std::vector<Space_Vector_>& cell_centers) {
	const auto num_cell = cell_centers.size();
	
	constexpr auto gamma = 1.4;
	constexpr auto c = 1 / (1.4 - 1);
	std::vector<Solution_> solutions(num_cell);
	for (size_t i = 0; i < num_cell; ++i) {
		const auto x_coordinate = cell_centers[i].at(0);

		if (x_coordinate <= 0.3) {
			constexpr auto rho = 1.0;
			constexpr auto u = 0.75;
			constexpr auto v = 0.0;
			constexpr auto p = 1.0;

			constexpr auto rhou = rho * u;
			constexpr auto rhov = rho * v;
			constexpr auto rhoE = p * c + 0.5 * (rhou * u + rhov * v);

			solutions[i] = { rho, rhou, rhov, rhoE };
		}
		else {
			constexpr auto rho = 0.125;
			constexpr auto u = 0.0;
			constexpr auto v = 0.0;
			constexpr auto p = 0.1;

			constexpr auto rhou = rho * u;
			constexpr auto rhov = rho * v;
			constexpr auto rhoE = p * c + 0.5 * (rhou * u + rhov * v);

			solutions[i] = { rho, rhou, rhov, rhoE };
		}
	}

	return solutions;
}

namespace ms {
	std::string d_to_str_with_precision(const double value, const ushort precision) {
		std::ostringstream os;
		os << std::setprecision(precision) << value;
		return os.str();
	}
}
