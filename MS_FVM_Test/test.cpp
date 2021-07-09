#include "gtest/gtest.h"

#include "../MS_FVM/RSC/Setting.h"
#include "../MS_FVM/INC/Governing_Equation.h"
#include "../MS_FVM/INC/Grid_Information_Builder.h"
#include "../MS_FVM/INC/Inital_Condition.h"

//#include <iomanip>

GTEST_TEST(TestBed, test){
	//const auto grid_data = GRID_READER::read("RSC/Grid/Quad_10.msh");
	//const auto grid_info = Grid_Information_Builder::build(grid_data);
	//const auto initial_solutions = INITIAL_CONDITION::calculate_initial_solutions(grid_info.cell_grid_information.centers);

	////for (size_t i = 0; i < 10; ++i) {
	////	for (size_t j = 0; j < 10; ++j) {
	////		std::cout << grid_info.cell_grid_information.centers[i * 10 + j] << "\t";
	////	}
	////	std::cout << "\n";
	////}


	//const auto num_solution = initial_solutions.size();
	//for (size_t i = 0; i < 10; ++i) {
	//	const auto y_coord = 0.05 + 0.1 * i;
	//	for (size_t j = 0; j < 10; ++j) {
	//		const auto x_coord = 0.05 + 0.1 * j;

	//		const auto result = grid_info.cell_grid_information.centers[i * 10 + j];

	//		Physical_Domain_Vector ref_ceter = { x_coord, y_coord };

	//		for(size_t i=0; i<PHYSICAL_DOMAIN_DIMENSION; ++i)
	//			EXPECT_DOUBLE_EQ(result[i], ref_ceter[i]);

	//		//ref_sol[i * 10 + j] = 2 * std::sin(2 * std::numbers::pi * x_coord) * std::sin(2 * std::numbers::pi * y_coord);
	//	}
	//}

	////EXPECT_EQ(grid_info.cell_grid_information.centers, ref_centers);
	////EXPECT_EQ(initial_solutions, ref_sol);
 }