#pragma once
#include "gtest/gtest.h"
#include "../MS_FVM/INC/EuclideanVector.h"

GTEST_TEST(EuclideanVector, operator_addition_1) {
	const EuclideanVector v1 = { 1,2,3 };
	const EuclideanVector v2 = { 4,5,6 };
	const auto result = v1 + v2;

	const EuclideanVector ref = { 5,7,9 };
	EXPECT_EQ(result, ref);
}
GTEST_TEST(EuclideanVector, operator_addition_2) {
	const EuclideanVector v1 = { 1,2,3,4 };
	const EuclideanVector v2 = { 4,5,6,7 };
	const auto result = v1 + v2;

	const EuclideanVector ref = { 5,7,9,11 };
	EXPECT_EQ(result, ref);
}

GTEST_TEST(EuclideanVector, operator_substraction_1) {
	const EuclideanVector v1 = { 1,2,3 };
	const EuclideanVector v2 = { 4,5,6 };
	const auto result = v1 - v2;

	const EuclideanVector ref = { -3,-3,-3 };
	EXPECT_EQ(result, ref);
}

GTEST_TEST(EuclideanVector, operator_scalar_multiplication_1) {
	const EuclideanVector v1 = { 1,2,3 };
	const auto result = v1 * 2;

	const EuclideanVector ref = { 2,4,6 };
	EXPECT_EQ(result, ref);
}
GTEST_TEST(EuclideanVector, operator_scalar_multiplication_2) {
	const EuclideanVector v1 = { 1,2,3 };
	const auto result = 2 * v1;

	const EuclideanVector ref = { 2,4,6 };
	EXPECT_EQ(result, ref);
}

GTEST_TEST(EuclideanVector, inner_product_1) {
	const EuclideanVector v1 = { 1,2,3 };
	const EuclideanVector v2 = { 4,5,6 };
	const auto result = v1.inner_product(v2);

	const double ref = 32;
	EXPECT_EQ(result, ref);
}

GTEST_TEST(EuclideanVector, norm_1) {
	const EuclideanVector v1 = { 3,4 };
	const auto result = v1.norm();

	const auto ref = 5;
	EXPECT_EQ(result, ref);
}