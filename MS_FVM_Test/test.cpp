#include "gtest/gtest.h"

#include "../MS_FVM/INC/Matrix.h"


TEST(TestCaseName, TestName) {
	Matrix<2, 2> m1;
	//Matrix<2, 2> m2 = { 1,2,3 };
	Matrix<2, 2> m3 = { 1,2,3,4 };

	std::vector<Matrix<2, 2>> mvec(6);
	for (size_t i = 0; i < 6; ++i) {
		mvec[i] = { 1,0,0,1 };
		//std::cout << mvec[i];
	}

}