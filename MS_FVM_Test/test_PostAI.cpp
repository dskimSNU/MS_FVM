//#pragma once
//#include "gtest/gtest.h"
//#include "../MS_FVM/INC/PostAI.h"
//
//using Grid_Builder_			= Grid_Builder<2>;
//
//
//GTEST_TEST(PostAI, calculate_face_share_cell_indexes_set1) {
//	auto grid = Grid_Builder_::build<Gmsh>("Quad3");
//
//	const auto face_share_cell_indexes_set = PostAI::calculate_face_share_cell_indexes_set(grid);
//
//	const auto result1 = face_share_cell_indexes_set[0];
//	const std::set<size_t> ref1 = { 1,2,3,6 };
//	EXPECT_EQ(result1, ref1);
//}
//GTEST_TEST(PostAI, calculate_face_share_cell_indexes_set2) {
//	auto grid = Grid_Builder_::build<Gmsh>("Quad3");
//
//	const auto face_share_cell_indexes_set = PostAI::calculate_face_share_cell_indexes_set(grid);
//
//	const auto result1 = face_share_cell_indexes_set[2];
//	const std::set<size_t> ref1 = { 0,1,5,8 };
//	EXPECT_EQ(result1, ref1);
//}
//GTEST_TEST(PostAI, calculate_face_share_cell_indexes_set3) {
//	auto grid = Grid_Builder_::build<Gmsh>("Quad3");
//
//	const auto face_share_cell_indexes_set = PostAI::calculate_face_share_cell_indexes_set(grid);
//
//	const auto result1 = face_share_cell_indexes_set[7];
//	const std::set<size_t> ref1 = { 1,4,6,8 };
//	EXPECT_EQ(result1, ref1);
//}
//
//GTEST_TEST(PostAI, calculate_vertex_nodes_coordinate_string_set) {
//	auto grid = Grid_Builder_::build<Gmsh>("Quad3");
//
//	const auto vnode_coorinate_string_set = PostAI::calculate_vertex_nodes_coordinate_string_set(grid);
//
//	std::string str;
//	for (const auto& s : vnode_coorinate_string_set)
//		str += s;
//	std::cout << str;
//}
//
//GTEST_TEST(PostAI, initialize) {
//	//auto grid_elements = Grid_File_Convertor_::convert_to_grid_elements("quad4");
//	//auto grid = Grid_Builder_::build(std::move(grid_elements));
//
//	//PostAI::intialize(grid);
//	//	
//	//std::cout << PostAI::ai_data_text_set_[5];
//
//	//const auto result = PostAI::ai_data_text_set_[5].size();
//	//const auto ref = 5;
//	//EXPECT_EQ(result, ref);
//}
//
//GTEST_TEST(PostAI, convert_to_solution_strings){
//	std::vector<Euclidean_Vector<1>> solutions;	
//	for (size_t i = 0; i < 10; ++i)
//		solutions.push_back({ i });
//	
//	const auto result = PostAI::convert_to_solution_strings(solutions);
//
//	for (const auto& str : result)
//		std::cout << str << "\n";
//}
//
//GTEST_TEST(PostAI, convert_to_solution_gradient_strings) {
//	std::vector<Dynamic_Matrix_> solution_gradients;
//	for (size_t i = 0; i < 10; ++i)
//		solution_gradients.push_back(Dynamic_Matrix_(1, 2, { 0.1 * i,0.2 * i }));
//
//	const auto result = PostAI::convert_to_solution_gradient_strings(solution_gradients);
//
//	for (const auto& str : result)
//		std::cout << str << "\n";
//}
//
////GTEST_TEST(PostAI, record_solution_data) {
////	auto grid_elements = Grid_File_Convertor_::convert_to_grid_elements("quad4");
////	auto grid = Grid_Builder_::build(std::move(grid_elements));
////
////	PostAI::intialize(grid);
////
////	constexpr auto num_solution = 16;
////	std::vector<Euclidean_Vector<1>> solutions;
////	for (size_t i = 0; i < num_solution; ++i)
////		solutions.push_back({ i });
////
////	std::vector<Dynamic_Matrix_> solution_gradients;
////	for (size_t i = 0; i < num_solution; ++i)
////		solution_gradients.push_back(Dynamic_Matrix_(1, 2, { 0.1 * i,-0.1 * i }));
////
////	PostAI::record_solution_datas(solutions, solution_gradients);
////
////	std::cout << PostAI::ai_data_text_set_[0];
////
////	const auto result = PostAI::ai_data_text_set_[0].size();
////	const auto ref = 7;
////	EXPECT_EQ(result, ref);
////}
//
//GTEST_TEST(PostAI, post) {
//	std::string path = "./AI_Data/";
//	PostAI::set_path(path);
//	
//	auto grid = Grid_Builder_::build<Gmsh>("Quad3");
//	PostAI::intialize(grid);
//
//	constexpr auto num_solution = 16;
//	std::vector<Euclidean_Vector<1>> solutions;
//	for (size_t i = 0; i < num_solution; ++i)
//		solutions.push_back({ i });
//
//	std::vector<Dynamic_Matrix_> solution_gradients;
//	for (size_t i = 0; i < num_solution; ++i)
//		solution_gradients.push_back(Dynamic_Matrix_(1, 2, { 0.1 * i,-0.1 * i }));
//
//	PostAI::record_solution_datas(solutions, solution_gradients);
//
//	std::array<double, 1> limiter_value = { 1 };
//	for (size_t i = 0; i < num_solution; ++i)
//		PostAI::record_limiting_value(i, limiter_value);
//
//	PostAI::post();
//}