#include "../RSC/Setting.h"
#include "../INC/Governing_Equation.h"
#include "../INC/Grid_Information_Builder.h"
#include "../INC/Inital_Condition.h"

int main(void) {
	const auto grid_data = GRID_READER::read(GRID_FILE_PATH);
	const auto grid_info = Grid_Information_Builder::build(grid_data);
	const auto initial_solutions = INITIAL_CONDITION::calculate_initial_solutions(grid_info.cell_grid_information.centers);
}