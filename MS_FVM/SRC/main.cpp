#include "../INC/Inital_Condition.h"
#include "../INC/Discrete_Equation.h"
#include "../INC/Setting.h"
#include "../INC/Post.h"
#include "../INC/Log.h"

using Post_						= Post<GOVERNING_EQUATION>;
using Grid_Builder_				= Grid_Builder<DIMENSION>;
using Semi_Discrete_Equation_	= Semi_Discrete_Equation<GOVERNING_EQUATION, SPATIAL_DISCRETE_METHOD, RECONSTRUCTION_METHOD, NUMERICAL_FLUX>;
using Discrete_Equation_		= Discrete_Equation<TIME_INTEGRAL_METHOD>;

int main(void) {
	Log::set_path(PATH);
	Post_::set_path(PATH);
	PostAI::set_path(PATH + "AI_Data/");

	auto grid = Grid_Builder_::build<GRID_FILE_TYPE>(GRID_FILE_NAME);

	Post_::intialize();
	Post_::grid(grid.elements.cell_elements);
	PostAI::intialize(grid);

	const auto semi_discrete_eq = Semi_Discrete_Equation_(std::move(grid));
	auto solutions				= semi_discrete_eq.calculate_initial_solutions<INITIAL_CONDITION>();
	
	Discrete_Equation_::solve<TIME_STEP_METHOD, SOLVE_END_CONDITION, SOLVE_POST_CONDITION, Post_>(semi_discrete_eq, solutions);
	semi_discrete_eq.estimate_error<INITIAL_CONDITION>(solutions, END_CONDITION_CONSTANT);

	Log::write();
}