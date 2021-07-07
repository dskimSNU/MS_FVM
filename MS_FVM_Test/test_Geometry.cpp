#pragma once
#include "gtest/gtest.h"
#include "../MS_FVM/INC/Geometry.h"

GTEST_TEST(ReferenceGeometry, quadrilateral_area_1){
	Physical_Domain_Vector n1 = { 1,1 };
	Physical_Domain_Vector n2 = { 1,2 };
	Physical_Domain_Vector n3 = { 2,2 };
	Physical_Domain_Vector n4 = { 2,1 };
	std::vector<Physical_Domain_Vector> nodes = { n1,n2,n3,n4 };

	ReferenceGeometry reference_geometry(ElementFigure::quadrilateral, 1);
	const auto result = reference_geometry.calculate_volume(nodes);

	const auto ref = 1;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(ReferenceGeometry, quadrilateral_area_2) {
	Physical_Domain_Vector n1 = { 1,1 };
	Physical_Domain_Vector n2 = { 2,1 };
	Physical_Domain_Vector n3 = { 4,2 };
	Physical_Domain_Vector n4 = { 1,2 };
	std::vector<Physical_Domain_Vector> nodes = { n1,n2,n3,n4 };

	ReferenceGeometry reference_geometry(ElementFigure::quadrilateral, 1);
	const auto result = reference_geometry.calculate_volume(nodes);

	const auto ref = 2;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(ReferenceGeometry, quadrilateral_area_3) {
	Physical_Domain_Vector n1 = { 1,1 };
	Physical_Domain_Vector n2 = { 2,1 };
	Physical_Domain_Vector n3 = { 4,2 };
	Physical_Domain_Vector n4 = { -100,1 };
	std::vector<Physical_Domain_Vector> nodes = { n1,n2,n3,n4 };

	ReferenceGeometry reference_geometry(ElementFigure::quadrilateral, 1);
	const auto result = reference_geometry.calculate_volume(nodes);

	const auto ref = 51;
	EXPECT_EQ(result, ref);
}

GTEST_TEST(ReferenceGeometry, triangle_area_1) {
	Physical_Domain_Vector n1 = { 1,1 };
	Physical_Domain_Vector n2 = { 2,1 };
	Physical_Domain_Vector n3 = { 4,2 };
	std::vector<Physical_Domain_Vector> nodes = { n1,n2,n3 };

	ReferenceGeometry reference_geometry(ElementFigure::triangle, 1);
	const auto result = reference_geometry.calculate_volume(nodes);

	const auto ref = 0.5;
	EXPECT_EQ(result, ref);	
}
GTEST_TEST(ReferenceGeometry, triangle_area_2) {
	Physical_Domain_Vector n1 = { 1.524,1 };
	Physical_Domain_Vector n2 = { 2,1.154 };
	Physical_Domain_Vector n3 = { 4.47411,2 };
	std::vector<Physical_Domain_Vector> nodes = { n1,n2,n3 };

	ReferenceGeometry reference_geometry(ElementFigure::triangle, 1);
	const auto result = reference_geometry.calculate_volume(nodes);

	const auto ref = 0.01084153;	
	//EXPECT_EQ(result, ref); //suffer by round off error
	//EXPECT_DOUBLE_EQ(result, ref); //suffer by round off error
	EXPECT_NEAR(result, ref, 1.0E-16);
}

GTEST_TEST(Geometry, faces_nodes_1) {
	ReferenceGeometry reference_geometry(ElementFigure::triangle, 1);

	Physical_Domain_Vector n1 = { 1.524,1 };
	Physical_Domain_Vector n2 = { 2,1.154 };
	Physical_Domain_Vector n3 = { 4.47411,2 };
	std::vector<Physical_Domain_Vector> nodes = { n1,n2,n3 };

	Geometry geometry(reference_geometry, std::move(nodes));
	const auto result = geometry.calculate_faces_nodes();

	const std::vector<std::vector<Physical_Domain_Vector>> ref = { {n1,n2},{n2,n3},{n3,n1} };
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Geometry, faces_nodes_2) {
	ReferenceGeometry reference_geometry(ElementFigure::quadrilateral, 1);

	Physical_Domain_Vector n1 = { 1,1 };
	Physical_Domain_Vector n2 = { 2,1 };
	Physical_Domain_Vector n3 = { 4,2 };
	Physical_Domain_Vector n4 = { 1,2 };
	std::vector<Physical_Domain_Vector> nodes = { n1,n2,n3,n4 };

	Geometry geometry(reference_geometry, std::move(nodes));
	const auto result = geometry.calculate_faces_nodes();

	const std::vector<std::vector<Physical_Domain_Vector>> ref = { {n1,n2},{n2,n3},{n3,n4}, {n4,n1} };
	EXPECT_EQ(result, ref);
}

GTEST_TEST(Geometry, coordinate_projected_volume_1) {
	ReferenceGeometry reference_geometry(ElementFigure::triangle, 1);

	Physical_Domain_Vector n1 = { 1.524,1 };
	Physical_Domain_Vector n2 = { 2,1.154 };
	Physical_Domain_Vector n3 = { 4.47411,2 };
	std::vector<Physical_Domain_Vector> nodes = { n1,n2,n3 };

	Geometry geometry(reference_geometry, std::move(nodes));
	const auto result = geometry.calculate_coordinate_projected_volume();

	const std::array<double, s_physical_domain_dimension> ref = { 4.47411 - 1.524, 1 };
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Geometry, coordinate_projected_volume_2) {
	ReferenceGeometry reference_geometry(ElementFigure::quadrilateral, 1);

	Physical_Domain_Vector n1 = { 1,1 };
	Physical_Domain_Vector n2 = { 2,1 };
	Physical_Domain_Vector n3 = { 4,2 };
	Physical_Domain_Vector n4 = { 1,2 };
	std::vector<Physical_Domain_Vector> nodes = { n1,n2,n3,n4 };

	Geometry geometry(reference_geometry, std::move(nodes));
	const auto result = geometry.calculate_coordinate_projected_volume();

	const std::array<double, s_physical_domain_dimension> ref = { 3,1 };
	EXPECT_EQ(result, ref);
}