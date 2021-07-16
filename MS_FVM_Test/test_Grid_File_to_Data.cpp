//#pragma once
//#include "gtest/gtest.h"
//#include "../MS_FVM/INC/Grid_File_to_Data.h"
//
//#include <random>
//
//GTEST_TEST(Grid_File_To_Data_Gmsh, read_about){
//	const std::string grid_file_path = "RSC/Grid/Mix_20.msh";
//	std::ifstream grid_file_stream(grid_file_path);
//	dynamic_require(grid_file_stream.is_open(), "fail to open grid file!");
//	
//	const auto node_text = Grid_File_Convertor<Gmsh, 2>::read_about(grid_file_stream, "Nodes");
//	//std::cout << node_text;
//
//	const auto element_text = Grid_File_Convertor<Gmsh, 2>::read_about(grid_file_stream, "Elements");
//	//std::cout << element_text;
//
//	const auto physical_name_text = Grid_File_Convertor<Gmsh, 2>::read_about(grid_file_stream, "PhysicalNames");
//	//std::cout << physical_name_text;	
//}
//
//GTEST_TEST(Grid_File_To_Data_Gmsh, make_node_data) {
//	const std::string grid_file_path = "RSC/Grid/Mix_20.msh";
//	std::ifstream grid_file_stream(grid_file_path);
//	dynamic_require(grid_file_stream.is_open(), "fail to open grid file!");
//
//	const auto node_text = Grid_File_Convertor<Gmsh, 2>::read_about(grid_file_stream, "Nodes");
//	const auto node_data_set = Grid_File_Convertor<Gmsh, 2>::make_node_grid_data(node_text);
//
//	for (size_t i = 0; i < node_data_set.size(); ++i) {
//		const auto result = node_data_set[i].to_string() + " ";
//		EXPECT_NE(node_text[i].find(result), std::string::npos);
//	}
//}
//
//GTEST_TEST(Grid_File_To_Data_Gmsh, physical_name) {
//	const std::string grid_file_path = "RSC/Grid/Mix_20.msh";
//	std::ifstream grid_file_stream(grid_file_path);
//	dynamic_require(grid_file_stream.is_open(), "fail to open grid file!");
//
//	const auto physical_name_text = Grid_File_Convertor<Gmsh, 2>::read_about(grid_file_stream, "PhysicalNames");
//	std::map<size_t, ElementType> result;
//	for (const auto& physical_name_sentence : physical_name_text) {
//		const char delimiter = ' ';
//		const auto parsed_sentence_set = ms::parse(physical_name_sentence, delimiter);
//
//		//const size_t dimension = parsed_sentence_set[0].toValue<size_t>();
//		const auto index = ms::string_to_value<size_t>(parsed_sentence_set[1]);
//		const auto name = ms::erase(parsed_sentence_set[2], "\"");
//		const auto element_type = ms::string_to_element_type(name);
//
//		result.emplace(index, element_type);
//	}
//
//	const std::map<size_t, ElementType> ref = { {1,ElementType::cell}, {2,ElementType::x_periodic}, {3,ElementType::y_periodic} };
//	EXPECT_EQ(result, ref);
//}
//
//bool operator==(const ElementGridData& ed1, const ElementGridData& ed2) {
//	return	ed1.figure == ed2.figure && ed1.figure_order == ed2.figure_order &&
//		ed1.index == ed2.index && ed1.node_indexes == ed2.node_indexes &&
//		ed1.type == ed2.type;
//}
//
//GTEST_TEST(Grid_File_To_Data_Gmsh, cell_data) {
//	const std::string grid_file_path = "RSC/Grid/Mix_20.msh";
//	std::ifstream grid_file_stream(grid_file_path);
//	dynamic_require(grid_file_stream.is_open(), "fail to open grid file!");
//
//	const auto element_text = Grid_File_Convertor<Gmsh, 2>::read_about(grid_file_stream, "Elements");
//	const auto physical_name_text = Grid_File_Convertor<Gmsh, 2>::read_about(grid_file_stream, "PhysicalNames");
//	const auto element_data_set = Grid_File_Convertor<Gmsh, 2>::make_element_data(element_text, physical_name_text);
//
//	const auto& cell_data_set = element_data_set[0];
//	const auto& cell_sample_result = cell_data_set.front();
//
//	const ElementGridData cell_ref_result = { 1, Figure::triangle, 1, ElementType::cell, {64, 367, 63} };
//	EXPECT_EQ(cell_sample_result, cell_ref_result);
//}
//
//GTEST_TEST(Grid_File_To_Data_Gmsh, boudnary_data) {
//	const std::string grid_file_path = "RSC/Grid/Mix_20.msh";
//	std::ifstream grid_file_stream(grid_file_path);
//	dynamic_require(grid_file_stream.is_open(), "fail to open grid file!");
//
//	const auto element_text = Grid_File_Convertor<Gmsh, 2>::read_about(grid_file_stream, "Elements");
//	const auto physical_name_text = Grid_File_Convertor<Gmsh, 2>::read_about(grid_file_stream, "PhysicalNames");
//	const auto element_data_set = Grid_File_Convertor<Gmsh, 2>::make_element_data(element_text, physical_name_text);
//
//	const auto& boundary_data_set = element_data_set[1];
//	const auto result = boundary_data_set.size();
//
//	const auto ref = 0;
//	EXPECT_EQ(result, ref);
//}
//
//GTEST_TEST(Grid_File_To_Data_Gmsh, periodic_data) {
//	const std::string grid_file_path = "RSC/Grid/Mix_20.msh";
//	std::ifstream grid_file_stream(grid_file_path);
//	dynamic_require(grid_file_stream.is_open(), "fail to open grid file!");
//
//	const auto element_text = Grid_File_Convertor<Gmsh, 2>::read_about(grid_file_stream, "Elements");
//	const auto physical_name_text = Grid_File_Convertor<Gmsh, 2>::read_about(grid_file_stream, "PhysicalNames");
//	const auto element_data_set = Grid_File_Convertor<Gmsh, 2>::make_element_data(element_text, physical_name_text);
//
//	const auto& periodic_boundary_data_set = element_data_set[2];	
//	const auto& periodic_boundary_sample_result = periodic_boundary_data_set.front();
//
//	const ElementGridData periodic_boundary_ref_result = { 651, Figure::line, 1, ElementType::x_periodic, {20,21} };
//	EXPECT_EQ(periodic_boundary_sample_result, periodic_boundary_ref_result);
//}