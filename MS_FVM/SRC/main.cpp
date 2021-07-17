//#include "../INC/Inital_Condition.h"
//#include "../INC/Discrete_Equation.h"
//#include "../INC/Setting.h"
//

#include "../INC/Setting.h"
#include "../INC/Grid_Builder.h"
#include "../INC/Post.h"
using Grid_File_Convertor_		= Grid_File_Convertor<GRID_FILE_TYPE, PHYSICAL_DOMAIN_DIMENSION>;
using Grid_Builder_				= Grid_Builder<PHYSICAL_DOMAIN_DIMENSION>;
////using Grid_Info_Extractor_		= Grid_Info_Extractor<SPATIAL_DISCRETE_METHOD, RECONSTRUCTION_METHOD, PHYSICAL_DOMAIN_DIMENSION>;
//using Semi_Discrete_Equation_	= Semi_Discrete_Equation<GOVERNING_EQUATION, SPATIAL_DISCRETE_METHOD, RECONSTRUCTION_METHOD, NUMERICAL_FLUX>;
//using Discrete_Equation_		= Discrete_Equation<TIME_INTEGRAL_METHOD>;
//
int main(void) {
//	//InitialCondtion<INITIAL_CONDITION_NAME,GOV_EQ>
//	Post::intialize<GOVERNING_EQUATION, INITIAL_CONDITION>(GRID_FILE_NAME);
//
	auto grid_elements	= Grid_File_Convertor_::convert_to_grid_elements(GRID_FILE_NAME);
	const auto grid		= Grid_Builder_::build(std::move(grid_elements));

	Post::grid(grid.elements.cell_elements);

//  const auto semi_dscrete_eq = Semi_Discrete_Equation_(processed_grid_data);
//	const Semi_Discrete_Equation_ semi_discrete_eq(processed_grid_data);
//	auto solutions = semi_discrete_eq.calculate_initial_solutions<INITIAL_CONDITION>();
//	
//	Discrete_Equation_::solve<TIME_STEP_METHOD, SOLVE_END_CONDITION, SOLVE_POST_CONDITION>(semi_discrete_eq, solutions);
//	semi_discrete_eq.estimate_error<INITIAL_CONDITION>(solutions, END_CONDITION_CONSTANT);
}