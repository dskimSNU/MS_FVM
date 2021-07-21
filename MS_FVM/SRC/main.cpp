#include "../INC/Inital_Condition.h"
#include "../INC/Discrete_Equation.h"
#include "../INC/Setting.h"
#include "../INC/Post.h"
#include "../INC/Log.h"

using Grid_File_Convertor_		= Grid_File_Convertor<GRID_FILE_TYPE, DIMENSION>;
using Grid_Builder_				= Grid_Builder<DIMENSION>;
using Semi_Discrete_Equation_	= Semi_Discrete_Equation<GOVERNING_EQUATION, SPATIAL_DISCRETE_METHOD, RECONSTRUCTION_METHOD, NUMERICAL_FLUX>;
using Discrete_Equation_		= Discrete_Equation<TIME_INTEGRAL_METHOD>;

int main(void) {
	const auto path = "E:/Code/Result/MS_FVM/" + GOVERNING_EQUATION::name() + "/" + INITIAL_CONDITION::name() + "/" + SPATIAL_DISCRETE_METHOD::name() + "_" + RECONSTRUCTION_METHOD::name()  + "/" + GRID_FILE_NAME + "/";
	Log::set_path(path);
	Post::set_path(path);

	auto grid_elements	= Grid_File_Convertor_::convert_to_grid_elements(GRID_FILE_NAME);
	auto grid			= Grid_Builder_::build(std::move(grid_elements));

	Post::intialize<GOVERNING_EQUATION>();
	Post::grid(grid.elements.cell_elements);

	const auto semi_discrete_eq = Semi_Discrete_Equation_(std::move(grid));
	auto solutions				= semi_discrete_eq.calculate_initial_solutions<INITIAL_CONDITION>();
	
	Discrete_Equation_::solve<TIME_STEP_METHOD, SOLVE_END_CONDITION, SOLVE_POST_CONDITION>(semi_discrete_eq, solutions);
	semi_discrete_eq.estimate_error<INITIAL_CONDITION>(solutions, END_CONDITION_CONSTANT);

	Log::write();
}