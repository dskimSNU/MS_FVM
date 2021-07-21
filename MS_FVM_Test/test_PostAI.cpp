#pragma once
#include "gtest/gtest.h"
#include "../MS_FVM/INC/Setting.h"
#include "../MS_FVM/INC/PostAI.h"

using Grid_File_Convertor_	= Grid_File_Convertor<GRID_FILE_TYPE, DIMENSION>;
using Grid_Builder_			= Grid_Builder<DIMENSION>;


GTEST_TEST(PostAI, calculate_face_share_cell_indexes_set1) {
	auto grid_elements = Grid_File_Convertor_::convert_to_grid_elements("quad3");
	auto grid = Grid_Builder_::build(std::move(grid_elements));

	const auto face_share_cell_indexes_set = PostAI::calculate_face_share_cell_indexes_set(grid);
	
	const auto result1 = face_share_cell_indexes_set[0];
	const std::set<size_t> ref1 = { 1,2,3,6 };
	EXPECT_EQ(result1, ref1);
}
GTEST_TEST(PostAI, calculate_face_share_cell_indexes_set2) {
	auto grid_elements = Grid_File_Convertor_::convert_to_grid_elements("quad3");
	auto grid = Grid_Builder_::build(std::move(grid_elements));

	const auto face_share_cell_indexes_set = PostAI::calculate_face_share_cell_indexes_set(grid);

	const auto result1 = face_share_cell_indexes_set[2];
	const std::set<size_t> ref1 = { 0,1,5,8 };
	EXPECT_EQ(result1, ref1);
}
GTEST_TEST(PostAI, calculate_face_share_cell_indexes_set3) {
	auto grid_elements = Grid_File_Convertor_::convert_to_grid_elements("quad3");
	auto grid = Grid_Builder_::build(std::move(grid_elements));

	const auto face_share_cell_indexes_set = PostAI::calculate_face_share_cell_indexes_set(grid);

	const auto result1 = face_share_cell_indexes_set[7];
	const std::set<size_t> ref1 = {1,4,6,8 };
	EXPECT_EQ(result1, ref1);
}

GTEST_TEST(PostAI, calculate_vertex_nodes_coordinate_string_set) {
	auto grid_elements = Grid_File_Convertor_::convert_to_grid_elements("quad3");
	auto grid = Grid_Builder_::build(std::move(grid_elements));

	const auto vnode_coorinate_string_set = PostAI::calculate_vertex_nodes_coordinate_string_set(grid);

	std::string str;
	for (const auto& s : vnode_coorinate_string_set)
		str += s;
	std::cout << str;
}

GTEST_TEST(PostAI, initialize)
{
	auto grid_elements = Grid_File_Convertor_::convert_to_grid_elements("quad4");
	auto grid = Grid_Builder_::build(std::move(grid_elements));

	PostAI::intialize(grid);

	std::cout << PostAI::node_number_string_set_[0];
	std::cout << PostAI::edge_number_string_set_[0];
	std::cout << PostAI::connectivity_string_set_[0];
	std::cout << PostAI::cell_coords_string_set_[0];
}