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
	std::vector<size_t> indexes = { 1,2,3 };

	Geometry geometry(reference_geometry, std::move(nodes), std::move(indexes));
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
	std::vector<size_t> indexes = { 1,2,3,4 };

	Geometry geometry(reference_geometry, std::move(nodes), std::move(indexes));
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
	std::vector<size_t> indexes = { 1,2,3 };

	Geometry geometry(reference_geometry, std::move(nodes), std::move(indexes));
	const auto result = geometry.coordinate_projected_volume();

	const std::array<double, Physical_Domain_Vector::dimension()> ref = { 4.47411 - 1.524, 1 };
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Geometry, coordinate_projected_volume_2) {
	ReferenceGeometry reference_geometry(ElementFigure::quadrilateral, 1);

	Physical_Domain_Vector n1 = { 1,1 };
	Physical_Domain_Vector n2 = { 2,1 };
	Physical_Domain_Vector n3 = { 4,2 };
	Physical_Domain_Vector n4 = { 1,2 };
	std::vector<Physical_Domain_Vector> nodes = { n1,n2,n3,n4 };
	std::vector<size_t> indexes = { 1,2,3,4 };

	Geometry geometry(reference_geometry, std::move(nodes), std::move(indexes));
	const auto result = geometry.coordinate_projected_volume();

	const std::array<double, Physical_Domain_Vector::dimension()> ref = { 3,1 };
	EXPECT_EQ(result, ref);
}

GTEST_TEST(Geometry, center_1) {
	ReferenceGeometry reference_geometry(ElementFigure::quadrilateral, 1);

	Physical_Domain_Vector n1 = { 1,1 };
	Physical_Domain_Vector n2 = { 2,1 };
	Physical_Domain_Vector n3 = { 4,2 };
	Physical_Domain_Vector n4 = { 1,2 };
	std::vector<Physical_Domain_Vector> nodes = { n1,n2,n3,n4 };
	std::vector<size_t> indexes = { 1,2,3,4 };

	Geometry geometry(reference_geometry, std::move(nodes), std::move(indexes));
	const auto result = geometry.center_node();

	const Physical_Domain_Vector ref = { 2,1.5 };
	EXPECT_EQ(result, ref);
}

GTEST_TEST(Geometry, calculate_normal_1) {
	ReferenceGeometry reference_geometry(ElementFigure::line, 1);

	Physical_Domain_Vector n1 = { 1,1 };
	Physical_Domain_Vector n2 = { 3,1 };
	std::vector<Physical_Domain_Vector> nodes = { n1,n2 };
	std::vector<size_t> indexes = { 1,2 };

	Geometry geometry(reference_geometry, std::move(nodes), std::move(indexes));

	Physical_Domain_Vector cell_center = { 2,0 };

	const auto normal = geometry.normal_vector(cell_center);

	const auto tangent = n2 - n1;
	const auto normality = tangent.inner_product(normal);
	const double size = normal.norm();
	const double direction = normal.inner_product({ 0,1 });


	const double normality_ref = 0;
	const double size_ref = 1;
	const double direction_ref = 1;

	EXPECT_EQ(normality, normality_ref);
	EXPECT_EQ(size, size_ref);
	EXPECT_EQ(direction, direction_ref);
}
GTEST_TEST(Geometry, calculate_normal_2) {
	ReferenceGeometry reference_geometry(ElementFigure::line, 1);

	Physical_Domain_Vector n1 = { 1,1 };
	Physical_Domain_Vector n2 = { 3,1 };
	std::vector<Physical_Domain_Vector> nodes = { n1,n2 };
	std::vector<size_t> indexes = { 1,2 };

	Geometry geometry(reference_geometry, std::move(nodes), std::move(indexes));

	Physical_Domain_Vector cell_center = { 2,2 };

	const auto normal = geometry.normal_vector(cell_center);

	const auto tangent = n2 - n1;
	const auto normality = tangent.inner_product(normal);
	const double size = normal.norm();
	const double direction = normal.inner_product({ 0,1 });
	
	//std::cout << normal; -0이 나오는데 뭔가 좀 이상하다~

	const double normality_ref = 0;
	const double size_ref = 1;
	const double direction_ref = -1;

	EXPECT_EQ(normality, normality_ref);
	EXPECT_EQ(size, size_ref);
	EXPECT_EQ(direction, direction_ref);
}

GTEST_TEST(Geometry, periodic_match_1) {
	ReferenceGeometry reference_geometry(ElementFigure::line, 1);

	Physical_Domain_Vector n1 = { 1,0 };
	Physical_Domain_Vector n2 = { 3,1 };
	std::vector<Physical_Domain_Vector> nodes = { n1,n2 };
	std::vector<size_t> indexes = { 1,2 };

	Geometry geometry(reference_geometry, std::move(nodes), std::move(indexes));

	Physical_Domain_Vector n3 = { 8,0 };
	Physical_Domain_Vector n4 = { 10,1 };
	std::vector<Physical_Domain_Vector> nodes2 = { n3,n4 };
	std::vector<size_t> indexes2 = { 1,2 };

	Geometry geometry2(reference_geometry, std::move(nodes2), std::move(indexes2));

	EXPECT_TRUE(geometry.is_periodic_pair(geometry2, ElementType::x_periodic));
}
GTEST_TEST(Geometry, periodic_match_2) {
	ReferenceGeometry reference_geometry(ElementFigure::line, 1);

	Physical_Domain_Vector n1 = { 1,0 };
	Physical_Domain_Vector n2 = { 3,1 };
	std::vector<Physical_Domain_Vector> nodes = { n1,n2 };
	std::vector<size_t> indexes = { 1,2 };

	Geometry geometry(reference_geometry, std::move(nodes), std::move(indexes));

	Physical_Domain_Vector n3 = { 8,0 };
	Physical_Domain_Vector n4 = { 10,1 };
	std::vector<Physical_Domain_Vector> nodes2 = { n3,n4 };
	std::vector<size_t> indexes2 = { 1,2 };

	Geometry geometry2(reference_geometry, std::move(nodes2), std::move(indexes2));

	EXPECT_FALSE(geometry.is_periodic_pair(geometry2, ElementType::y_periodic));
}
GTEST_TEST(Geometry, periodic_match_3) {
	ReferenceGeometry reference_geometry(ElementFigure::line, 1);

	Physical_Domain_Vector n1 = { 1,0 };
	Physical_Domain_Vector n2 = { 3,1 };
	std::vector<Physical_Domain_Vector> nodes = { n1,n2 };
	std::vector<size_t> indexes = { 1,2 };

	Geometry geometry(reference_geometry, std::move(nodes), std::move(indexes));

	Physical_Domain_Vector n3 = { 1,5 };
	Physical_Domain_Vector n4 = { 3,6 };
	std::vector<Physical_Domain_Vector> nodes2 = { n3,n4 };
	std::vector<size_t> indexes2 = { 1,2 };

	Geometry geometry2(reference_geometry, std::move(nodes2), std::move(indexes2));

	EXPECT_TRUE(geometry.is_periodic_pair(geometry2, ElementType::y_periodic));
}

GTEST_TEST(Geometry, faces_geometry_1) {
	ReferenceGeometry reference_geometry(ElementFigure::triangle, 1);

	Physical_Domain_Vector n1 = { 1,1 };
	Physical_Domain_Vector n2 = { 2,1 };
	Physical_Domain_Vector n3 = { 4,2 };
	std::vector<Physical_Domain_Vector> nodes = { n1,n2,n3 };
	std::vector<size_t> indexes = { 1,2,3 };

	Geometry geometry(reference_geometry, std::move(nodes), std::move(indexes));
	const auto result = geometry.face_geometries();

	ReferenceGeometry reference_face_geometry(ElementFigure::line, 1);
	Geometry f1(reference_face_geometry, { {1,1}, {2,1} }, { 1,2 });
	Geometry f2(reference_face_geometry, { {2,1}, {4,2} }, { 1,2 });
	Geometry f3(reference_face_geometry, { {4,2}, {1,1} }, { 1,2 });

	const std::vector<Geometry> ref = { f1,f2,f3 };
	EXPECT_EQ(result, ref);
}