#include "../INC/Inital_Condition.h"
#include "../INC/Discrete_Equation.h"
#include "../INC/Setting.h"

using Grid_File_Reader_			= Grid_File_Reader<GRID_FILE_TYPE, PHYSICAL_DOMAIN_DIMENSION>;
using Grid_Data_Processor_		= Grid_Data_Processor<PHYSICAL_DOMAIN_DIMENSION>;
using Grid_Info_Extractor_		= Grid_Info_Extractor<SPATIAL_DISCRETE_METHOD, PHYSICAL_DOMAIN_DIMENSION>;
using Semi_Discrete_Equation_	= Semi_Discrete_Equation<SPATIAL_DISCRETE_METHOD, TIME_STEP_METHOD, GOVERNING_EQUATION, NUMERICAL_FLUX>;
using Discrete_Equation_		= Discrete_Equation<TIME_INTGRAL_METHOD>;

int main(void) {
	//InitialCondtion<INITIAL_CONDITION_NAME,GOV_EQ>
	Post::intialize<GOVERNING_EQUATION, INITIAL_CONDITION>(GRID_FILE_NAME);

	auto grid_raw_data			= Grid_File_Reader_::read(GRID_FILE_NAME);
	auto processed_grid_data	= Grid_Data_Processor_::process(std::move(grid_raw_data));
	auto grid_info				= Grid_Info_Extractor_::extract(std::move(processed_grid_data));
	const Semi_Discrete_Equation_ semi_discrete_eq(std::move(grid_info));

	auto solutions				= INITIAL_CONDITION::calculate_solutions(grid_info.cell_informations.centers);
	Discrete_Equation_::solve<SOLVE_END_CONDITION, SOLVE_POST_CONDITION>(semi_discrete_eq, solutions);
	semi_discrete_eq.estimate_error<INITIAL_CONDITION>(solutions, END_CONDITION_CONSTANT);
}