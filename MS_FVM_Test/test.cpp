#include "gtest/gtest.h"

#include "../MS_FVM/INC/Matrix.h"
#include "../MS_FVM/INC/Text.h"
#include "../MS_FVM/INC/Grid_Reader.h"

//#include <iomanip>

GTEST_TEST(TestBed, test){
	std::string str = "0.9443476362451639";

	std::istringstream iss(str);
	double value1;
	iss >> value1;

	const auto value2 = std::stod(str);
	EXPECT_EQ(value1, value2);
}