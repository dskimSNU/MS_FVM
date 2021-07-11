#include "gtest/gtest.h"
#include "../MS_FVM/INC/Grid_File_Type.h"

GTEST_TEST(Gmsh, figure_type_index_to_element_figure_1) {
	constexpr size_t figure_type_index = 1;
	const auto result = Gmsh::figure_type_index_to_element_figure(figure_type_index);

	const auto ref = Figure::line;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Gmsh, figure_type_index_to_element_figure_2) {
	constexpr size_t figure_type_index = 2;
	const auto result = Gmsh::figure_type_index_to_element_figure(figure_type_index);

	const auto ref = Figure::triangle;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Gmsh, figure_type_index_to_element_figure_3) {
	constexpr size_t figure_type_index = 3;
	const auto result = Gmsh::figure_type_index_to_element_figure(figure_type_index);

	const auto ref = Figure::quadrilateral;
	EXPECT_EQ(result, ref);
}

GTEST_TEST(Gmsh, figure_type_index_to_figure_order_1) {
	constexpr size_t figure_type_index = 1;
	const auto result = Gmsh::figure_type_index_to_figure_order(figure_type_index);

	const auto ref = 1;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Gmsh, figure_type_index_to_figure_order_2) {
	constexpr size_t figure_type_index = 2;
	const auto result = Gmsh::figure_type_index_to_figure_order(figure_type_index);

	const auto ref = 1;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Gmsh, figure_type_index_to_figure_order_3) {
	constexpr size_t figure_type_index = 3;
	const auto result = Gmsh::figure_type_index_to_figure_order(figure_type_index);

	const auto ref = 1;
	EXPECT_EQ(result, ref);
}