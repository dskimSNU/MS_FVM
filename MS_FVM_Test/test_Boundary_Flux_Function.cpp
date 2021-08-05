#pragma once
#include "gtest/gtest.h"
#include "../MS_FVM/INC/Boundary_Flux_Function.h"

GTEST_TEST(Supersonic_Outlet_2D, Linear_Advection_2D_1) {
	const auto boundary_type = ElementType::supersonic_outlet_2D;
	const auto bff = Boundary_Flux_Function_Factory<Linear_Advection_2D>::make(boundary_type);				

	const EuclideanVector solution = { 2 };
	const EuclideanVector normal = { 1,0 };
	const auto result = bff->calculate(solution, normal);

	const EuclideanVector ref = { 2 };
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Supersonic_Outlet_2D, Linear_Advection_2D_2) {
	const auto boundary_type = ElementType::supersonic_outlet_2D;
	const auto bff = Boundary_Flux_Function_Factory<Linear_Advection_2D>::make(boundary_type);

	const EuclideanVector solution = { 2 };
	const EuclideanVector normal = { 1,1 };
	const auto result = bff->calculate(solution, normal);

	const EuclideanVector ref = { 3 };
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Supersonic_Outlet_2D, Linear_Advection_2D_3) {
	const auto boundary_type = ElementType::supersonic_outlet_2D;
	const auto bff = Boundary_Flux_Function_Factory<Linear_Advection_2D>::make(boundary_type);

	const EuclideanVector solution = { 2 };
	const EuclideanVector normal = { 1,-1 };
	const auto result = bff->calculate(solution, normal);

	const EuclideanVector ref = { 1 };
	EXPECT_EQ(result, ref);
}


GTEST_TEST(Supersonic_Outlet_2D, Burgers_2D_1) {
	const auto boundary_type = ElementType::supersonic_outlet_2D;
	const auto bff = Boundary_Flux_Function_Factory<Burgers_2D>::make(boundary_type);

	const EuclideanVector solution = { 4 };
	const EuclideanVector normal = { 1,0 };
	const auto result = bff->calculate(solution, normal);

	const EuclideanVector ref = { 8 };
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Supersonic_Outlet_2D, Burgers_2D_2) {
	const auto boundary_type = ElementType::supersonic_outlet_2D;
	const auto bff = Boundary_Flux_Function_Factory<Burgers_2D>::make(boundary_type);

	const EuclideanVector solution = { 4 };
	const EuclideanVector normal = { 1,1 };
	const auto result = bff->calculate(solution, normal);

	const EuclideanVector ref = { 16 };
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Supersonic_Outlet_2D, Burgers_2D_3) {
	const auto boundary_type = ElementType::supersonic_outlet_2D;
	const auto bff = Boundary_Flux_Function_Factory<Burgers_2D>::make(boundary_type);

	const EuclideanVector solution = { 4 };
	const EuclideanVector normal = { 1,-1 };
	const auto result = bff->calculate(solution, normal);

	const EuclideanVector ref = { 0 };
	EXPECT_EQ(result, ref);
}


GTEST_TEST(Supersonic_Outlet_2D, Euler_2D_1) {
	const auto boundary_type = ElementType::supersonic_outlet_2D;
	const auto bff = Boundary_Flux_Function_Factory<Euler_2D>::make(boundary_type);

	const EuclideanVector cvariable = { 1,0,0,1 };
	const EuclideanVector normal = { 1,0 };
	const auto result = bff->calculate(cvariable, normal);

	const EuclideanVector ref = { 0, (1.4 - 1), 0, 0 };
	EXPECT_EQ(result, ref);
}
