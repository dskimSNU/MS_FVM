#pragma once
#include "gtest/gtest.h"
#include "../MS_FVM/INC/Governing_Equation.h"

GTEST_TEST(Linear_Advection_2D, calculate_physical_fluxes_1) {
	constexpr size_t num = 5;
	
	std::vector<Linear_Advection_2D::Solution> solutions(num);
	for (size_t i = 0; i < num; ++i)
		solutions[i] = { i };	
	const auto result = Linear_Advection_2D::calculate_physical_fluxes(solutions);

	std::vector<Matrix<1, 2>> ref(num);
	for (size_t i = 0; i < num; ++i)
		ref[i] = { i, 0.5 * i };

	EXPECT_EQ(result, ref);
}

GTEST_TEST(Linear_Advection_2D, calculate_coordinate_projected_maximum_lambdas_1) {
	constexpr size_t num = 5;

	std::vector<Linear_Advection_2D::Solution> solutions(num);
	for (size_t i = 0; i < num; ++i)
		solutions[i] = { i };
	const auto result = Linear_Advection_2D::calculate_coordinate_projected_maximum_lambdas(solutions);

	std::vector<double> x_projected_maximum_lambdas(num, 1.0);
	std::vector<double> y_projected_maximum_lambdas(num, 0.5);
	std::array<std::vector<double>, 2> ref = { x_projected_maximum_lambdas,y_projected_maximum_lambdas };
	EXPECT_EQ(result, ref);
}

GTEST_TEST(Linear_Advection_2D, calculate_inner_face_maximum_lambdas_1) {	
	Linear_Advection_2D::Solution solution_o = 1;
	Linear_Advection_2D::Solution solution_n = 0;

	constexpr size_t num = 5;
	for (size_t i = 0; i < num; ++i)
		Linear_Advection_2D::Physical_Domain_Vector normal = { i, i };

	std::vector<Linear_Advection_2D::Solution> solutions(num);
	for (size_t i = 0; i < num; ++i)
		solutions[i] = { i };
	const auto result = Linear_Advection_2D::calculate_coordinate_projected_maximum_lambdas(solutions);

	std::vector<double> x_projected_maximum_lambdas(num, 1.0);
	std::vector<double> y_projected_maximum_lambdas(num, 0.5);
	std::array<std::vector<double>, 2> ref = { x_projected_maximum_lambdas,y_projected_maximum_lambdas };
	EXPECT_EQ(result, ref);
}