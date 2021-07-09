#pragma once
#include "gtest/gtest.h"

#include "../MS_FVM/INC/Grid_Information_Builder.h"


GTEST_TEST(Grid_Data_Converter, cell_center) {
	const auto grid_data = Gmsh_Grid_Reader::read("RSC/Grid/Quad_10.msh");
	const auto grid_info = Grid_Data_Converter::convert(grid_data);

	const auto& cell_centers = grid_info.cell_grid_information.centers;

	const auto num_cell = cell_centers.size();
	for (size_t i = 0; i < 10; ++i) {
		const auto y_coord = 0.05 + 0.1 * i;
		for (size_t j = 0; j < 10; ++j) {
			const auto x_coord = 0.05 + 0.1 * j;

			const auto result = cell_centers[i * 10 + j];

			Physical_Domain_Vector ref_ceter = { x_coord, y_coord };
			for (size_t i = 0; i < PHYSICAL_DOMAIN_DIMENSION; ++i)
				EXPECT_DOUBLE_EQ(result[i], ref_ceter[i]);
		}
	}
}

GTEST_TEST(Grid_Data_Converter, volume) {
	const auto grid_data = Gmsh_Grid_Reader::read("RSC/Grid/Quad_10.msh");
	const auto grid_info = Grid_Data_Converter::convert(grid_data);

	const auto& cell_volumes = grid_info.cell_grid_information.volumes;
	
	const auto ref = 0.01;
	for (const auto& volume : cell_volumes)
		EXPECT_NEAR(volume, ref, 1.0E-16);
		//EXPECT_DOUBLE_EQ(volume, ref);	//suffer by round off error
}

GTEST_TEST(Grid_Data_Converter, coordinate_projected_volumes) {
	const auto grid_data = Gmsh_Grid_Reader::read("RSC/Grid/Quad_10.msh");
	const auto grid_info = Grid_Data_Converter::convert(grid_data);

	const auto& coordinate_projected_volumes = grid_info.cell_grid_information.coordinate_projected_volumes;

	const auto ref = 0.1;
	for (const auto& coordinate_projected_volumes : coordinate_projected_volumes) {
		const auto [x_projected, y_projected] = coordinate_projected_volumes;
		EXPECT_NEAR(x_projected, ref, 9.0E-16);
		EXPECT_NEAR(y_projected, ref, 9.0E-16);	//suffer by round off error
	}
}