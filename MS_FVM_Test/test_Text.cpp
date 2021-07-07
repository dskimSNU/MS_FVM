#include "gtest/gtest.h"
#include "../MS_FVM/INC/Text.h"


GTEST_TEST(ms, upper_case_1) {
	std::string str = "abc";
	const auto result = ms::upper_case(str);

	std::string ref = "ABC";
	EXPECT_EQ(result, ref);
}
GTEST_TEST(ms, upper_case_2) {
	std::string str = "abc123";
	const auto result = ms::upper_case(str);

	std::string ref = "ABC123";
	EXPECT_EQ(result, ref);
}
GTEST_TEST(ms, upper_case_3) {
	std::string str = "abc_123q";
	const auto result = ms::upper_case(str);

	std::string ref = "ABC_123Q";
	EXPECT_EQ(result, ref);
}

GTEST_TEST(ms, find_icase_1) {
	std::string str = "abc_123q";
	const auto result = ms::find_icase(str, "BC_");

	const auto ref = 1;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(ms, find_icase_2) {
	std::string str = "abc_123q";
	const auto result = ms::find_icase(str, "C_12");

	const auto ref = 2;
	EXPECT_EQ(result, ref);
}

GTEST_TEST(ms, is_there_icase_1) {
	std::string str = "abc_123q";
	const auto result = ms::is_there_icase(str, "C_12");

	const auto ref = true;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(ms, is_there_icase_2) {
	std::string str = "abc_123q";
	const auto result = ms::is_there_icase(str, "3Q");

	const auto ref = true;
	EXPECT_EQ(result, ref);
}