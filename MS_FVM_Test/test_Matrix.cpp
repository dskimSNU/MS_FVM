#pragma once
#include "gtest/gtest.h"
#include "../MS_FVM/INC/Matrix.h"

GTEST_TEST(Matrix, operator_mv_1) {
	const EuclideanVector v = { 1,2 };
	const Matrix<2,2> m = { 1,2,3,4 };
	const auto result = m * v;

	const EuclideanVector ref = { 5,11 };
	EXPECT_EQ(result, ref);
}