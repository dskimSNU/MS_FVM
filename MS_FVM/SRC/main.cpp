#include "../RSC/Setting.h"
#include "../INC/Governing_Equation.h"
#include "../INC/Grid_Information_Builder.h"
#include "../INC/Inital_Condition.h"
#include "../INC/Discrete_Equation.h"
#include "../INC/Time_Integral_Method.h"
#include "../INC/Solve_Condition.h"

int main(void) {
	const auto grid_data = GRID_READER::read(GRID_FILE_PATH);
	auto grid_info = Grid_Data_Converter::convert(grid_data);
	auto initial_solutions = INITIAL_CONDITION::calculate_solutions(grid_info.cell_grid_information.centers);

	SEMI_DISCRETE_EQUATION<GOVERNING_EQUATION, NUMERICAL_FLUX> semi_discrete_Eq(std::move(grid_info));
	Discrete_Equation::solve<TIME_INTGRAL_METHOD, END_CONDITION, POST_CONDITION>(semi_discrete_Eq, initial_solutions);

	//cell inner faces
	//semi discrete equation << cell inner face
	//discrete_eq << semi discrete equation <Time integral Method>;
	//solve<end,post condition>
}