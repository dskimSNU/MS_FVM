#include "gtest/gtest.h"

#include "../MS_FVM/RSC/Setting.h"
#include "../MS_FVM/INC/Governing_Equation.h"
#include "../MS_FVM/INC/Grid_Information_Builder.h"
#include "../MS_FVM/INC/Inital_Condition.h"
#include "../MS_FVM/INC/Discrete_Equation.h"
#include "../MS_FVM/INC/Time_Integral_Method.h"
#include "../MS_FVM/INC/Solve_Condition.h"

//#include <iomanip>

GTEST_TEST(TestBed, test){
	const auto grid_data = GRID_READER::read(GRID_FILE_PATH);
	auto grid_info = Grid_Data_Converter::convert(grid_data);
	auto initial_solutions = INITIAL_CONDITION::calculate_solutions(grid_info.cell_grid_information.centers);

	SEMI_DISCRETE_EQUATION<GOVERNING_EQUATION, NUMERICAL_FLUX> semi_discrete_Eq(std::move(grid_info));
	Discrete_Equation::solve<TIME_INTGRAL_METHOD, END_CONDITION, POST_CONDITION>(semi_discrete_Eq, initial_solutions);
 }