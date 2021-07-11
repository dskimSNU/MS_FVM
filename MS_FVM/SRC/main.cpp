#include "../RSC/Setting.h"
#include "../INC/Inital_Condition.h"
#include "../INC/Discrete_Equation.h"
#include "../INC/Post.h"

int main(void) {
	Post::intialize<GOVERNING_EQUATION, PHYSICAL_DOMAIN_DIMENSION>(POST_FILE_NAME);

	auto grid_data = Grid_File_To_Data<GRID_FILE_TYPE, PHYSICAL_DOMAIN_DIMENSION>::convert(GRID_FILE_PATH);
	auto grid_info = Grid_Data_to_Info<PHYSICAL_DOMAIN_DIMENSION>::convert(std::move(grid_data));
	auto initial_solutions = INITIAL_CONDITION::calculate_solutions(grid_info.cell_grid_information.centers);

	Semi_Discrete_Equation<SPATIAL_DISCRETE_METHOD, GOVERNING_EQUATION, NUMERICAL_FLUX> semi_discrete_Eq(std::move(grid_info));
	Discrete_Equation<TIME_INTGRAL_METHOD>::solve<END_CONDITION, POST_CONDITION>(semi_discrete_Eq, initial_solutions);
}