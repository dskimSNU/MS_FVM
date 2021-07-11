#pragma once
#include "gtest/gtest.h"
#include "../MS_FVM/INC/Geometry.h"

GTEST_TEST(Geometry, volume_1){
	const Figure fig = Figure::quadrilateral;
	const size_t fig_order = 1;

	const EuclideanVector n1 = { 1,1 };
	const EuclideanVector n2 = { 1,2 };
	const EuclideanVector n3 = { 2,2 };
	const EuclideanVector n4 = { 2,1 };
	std::vector<EuclideanVector<2>> nodes = { n1,n2,n3,n4 };
	std::vector<size_t> indexes = { 1,2,3,4 };

	Geometry geometry(fig, fig_order, std::move(nodes), std::move(indexes));
	const auto result = geometry.volume();

	const auto ref = 1;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Geometry, volume_2) {
	const Figure fig = Figure::quadrilateral;
	const size_t fig_order = 1;

	const EuclideanVector n1 = { 1,1 };
	const EuclideanVector n2 = { 2,1 };
	const EuclideanVector n3 = { 4,2 };
	const EuclideanVector n4 = { 1,2 };
	std::vector<EuclideanVector<2>> nodes = { n1,n2,n3,n4 };
	std::vector<size_t> indexes = { 1,2,3,4 };

	Geometry geometry(fig, fig_order, std::move(nodes), std::move(indexes));
	const auto result = geometry.volume();

	const auto ref = 2;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(ReferenceGeometry, volume_3) {
	const Figure fig = Figure::quadrilateral;
	const size_t fig_order = 1;

	const EuclideanVector n1 = { 1,1 };
	const EuclideanVector n2 = { 2,1 };
	const EuclideanVector n3 = { 4,2 };
	const EuclideanVector n4 = { -100,1 };
	std::vector<EuclideanVector<2>> nodes = { n1,n2,n3,n4 };
	std::vector<size_t> indexes = { 1,2,3,4 };

	Geometry geometry(fig, fig_order, std::move(nodes), std::move(indexes));
	const auto result = geometry.volume();

	const auto ref = 51;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(ReferenceGeometry, volume_4) {
	const Figure fig = Figure::triangle;
	const size_t fig_order = 1;

	const EuclideanVector n1 = { 1,1 };
	const EuclideanVector n2 = { 2,1 };
	const EuclideanVector n3 = { 4,2 };
	std::vector<EuclideanVector<2>> nodes = { n1,n2,n3 };
	std::vector<size_t> indexes = { 1,2,3,4 };

	Geometry geometry(fig, fig_order, std::move(nodes), std::move(indexes));
	const auto result = geometry.volume();

	const auto ref = 0.5;
	EXPECT_EQ(result, ref);	
}
GTEST_TEST(ReferenceGeometry, volume_5) {
	const Figure fig = Figure::triangle;
	const size_t fig_order = 1;

	const EuclideanVector n1 = { 1.524,1 };
	const EuclideanVector n2 = { 2,1.154 };
	const EuclideanVector n3 = { 4.47411,2 };

	std::vector<EuclideanVector<2>> nodes = { n1,n2,n3 };
	std::vector<size_t> indexes = { 1,2,3,4 };

	Geometry geometry(fig, fig_order, std::move(nodes), std::move(indexes));
	const auto result = geometry.volume();

	const auto ref = 0.01084153;	
	//EXPECT_EQ(result, ref); //suffer by round off error
	//EXPECT_DOUBLE_EQ(result, ref); //suffer by round off error
	EXPECT_NEAR(result, ref, 1.0E-16);
}

//GTEST_TEST(Geometry, faces_nodes_1) {
//	const Figure fig = Figure::triangle;
//	const size_t fig_order = 1;
//
//	const EuclideanVector n1 = { 1.524,1 };
//	const EuclideanVector n2 = { 2,1.154 };
//	const EuclideanVector n3 = { 4.47411,2 };
//
//	std::vector<EuclideanVector<2>> nodes = { n1,n2,n3 };
//	std::vector<size_t> indexes = { 1,2,3,4 };
//
//	Geometry geometry(fig, fig_order, std::move(nodes), std::move(indexes));
//	const auto result = geometry.calculate_faces_nodes();
//
//	const std::vector<std::vector<EuclideanVector>> ref = { {n1,n2},{n2,n3},{n3,n1} };
//	EXPECT_EQ(result, ref);
//}
//GTEST_TEST(Geometry, faces_nodes_2) {
//	ReferenceGeometry reference_geometry(Figure::quadrilateral, 1);
//
//	EuclideanVector n1 = { 1,1 };
//	EuclideanVector n2 = { 2,1 };
//	EuclideanVector n3 = { 4,2 };
//	EuclideanVector n4 = { 1,2 };
//	std::vector<EuclideanVector> nodes = { n1,n2,n3,n4 };
//	std::vector<size_t> indexes = { 1,2,3,4 };
//
//	Geometry geometry(reference_geometry, std::move(nodes), std::move(indexes));
//	const auto result = geometry.calculate_faces_nodes();
//
//	const std::vector<std::vector<EuclideanVector>> ref = { {n1,n2},{n2,n3},{n3,n4}, {n4,n1} };
//	EXPECT_EQ(result, ref);
//}

GTEST_TEST(Geometry, coordinate_projected_volume_1) {
	const Figure fig = Figure::triangle;
	const size_t fig_order = 1;

	const EuclideanVector n1 = { 1.524,1 };
	const EuclideanVector n2 = { 2,1.154 };
	const EuclideanVector n3 = { 4.47411,2 };

	std::vector<EuclideanVector<2>> nodes = { n1,n2,n3 };
	std::vector<size_t> indexes = { 1,2,3,4 };

	Geometry geometry(fig, fig_order, std::move(nodes), std::move(indexes));
	const auto result = geometry.coordinate_projected_volume();

	const std::array<double, 2> ref = { 4.47411 - 1.524, 1 };
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Geometry, coordinate_projected_volume_2) {
	const Figure fig = Figure::quadrilateral;
	const size_t fig_order = 1;

	const EuclideanVector n1 = { 1,1 };
	const EuclideanVector n2 = { 2,1 };
	const EuclideanVector n3 = { 4,2 };
	const EuclideanVector n4 = { 1,2 };
	std::vector<EuclideanVector<2>> nodes = { n1,n2,n3,n4 };
	std::vector<size_t> indexes = { 1,2,3,4 };

	Geometry geometry(fig, fig_order, std::move(nodes), std::move(indexes));
	const auto result = geometry.coordinate_projected_volume();

	const std::array<double, 2> ref = { 3,1 };
	EXPECT_EQ(result, ref);
}

GTEST_TEST(Geometry, center_1) {
	const Figure fig = Figure::quadrilateral;
	const size_t fig_order = 1;

	const EuclideanVector n1 = { 1,1 };
	const EuclideanVector n2 = { 2,1 };
	const EuclideanVector n3 = { 4,2 };
	const EuclideanVector n4 = { 1,2 };
	std::vector<EuclideanVector<2>> nodes = { n1,n2,n3,n4 };
	std::vector<size_t> indexes = { 1,2,3,4 };

	Geometry geometry(fig, fig_order, std::move(nodes), std::move(indexes));
	const auto result = geometry.center_node();

	const EuclideanVector ref = { 2,1.5 };
	EXPECT_EQ(result, ref);
}

GTEST_TEST(Geometry, calculate_normal_1) {
	const Figure fig = Figure::line;
	const size_t fig_order = 1;

	const EuclideanVector n1 = { 1,1 };
	const EuclideanVector n2 = { 3,1 };
	std::vector<EuclideanVector<2>> nodes = { n1,n2 };
	std::vector<size_t> indexes = { 1,2 };

	Geometry geometry(fig, fig_order, std::move(nodes), std::move(indexes));

	const EuclideanVector cell_center = { 2,0 };

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
	const Figure fig = Figure::line;
	const size_t fig_order = 1;

	const EuclideanVector n1 = { 1,1 };
	const EuclideanVector n2 = { 3,1 };
	std::vector<EuclideanVector<2>> nodes = { n1,n2 };
	std::vector<size_t> indexes = { 1,2 };

	Geometry geometry(fig, fig_order, std::move(nodes), std::move(indexes));

	EuclideanVector cell_center = { 2,2 };

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
	const Figure fig = Figure::line;
	const size_t fig_order = 1;

	const EuclideanVector n1 = { 1,1 };
	const EuclideanVector n2 = { 3,1 };
	std::vector<EuclideanVector<2>> nodes = { n1,n2 };
	std::vector<size_t> indexes = { 1,2 };

	Geometry geometry(fig, fig_order, std::move(nodes), std::move(indexes));
	
	const EuclideanVector n3 = { 8,1 };
	const EuclideanVector n4 = { 10,1 };
	std::vector<EuclideanVector<2>> nodes2 = { n3,n4 };
	std::vector<size_t> indexes2 = { 1,2 };

	Geometry geometry2(fig, fig_order, std::move(nodes2), std::move(indexes2));

	constexpr size_t axis_tag = 0;
	EXPECT_TRUE(geometry.is_periodic_pair(geometry2, axis_tag));
}
GTEST_TEST(Geometry, periodic_match_2) {
	const Figure fig = Figure::line;
	const size_t fig_order = 1;

	const EuclideanVector n1 = { 1,0 };
	const EuclideanVector n2 = { 3,1 };
	std::vector<EuclideanVector<2>> nodes = { n1,n2 };
	std::vector<size_t> indexes = { 1,2 };

	Geometry geometry(fig, fig_order, std::move(nodes), std::move(indexes));

	const EuclideanVector n3 = { 8,0 };
	const EuclideanVector n4 = { 10,1 };
	std::vector<EuclideanVector<2>> nodes2 = { n3,n4 };
	std::vector<size_t> indexes2 = { 1,2 };

	Geometry geometry2(fig, fig_order, std::move(nodes2), std::move(indexes2));

	constexpr size_t axis_tag = 1;
	EXPECT_FALSE(geometry.is_periodic_pair(geometry2, axis_tag));
}
GTEST_TEST(Geometry, periodic_match_3) {
	const Figure fig = Figure::line;
	const size_t fig_order = 1;

	const EuclideanVector n1 = { 1,0 };
	const EuclideanVector n2 = { 3,1 };
	std::vector<EuclideanVector<2>> nodes = { n1,n2 };
	std::vector<size_t> indexes = { 1,2 };

	Geometry geometry(fig, fig_order, std::move(nodes), std::move(indexes));

	const EuclideanVector n3 = { 1,5 };
	const EuclideanVector n4 = { 3,6 };
	std::vector<EuclideanVector<2>> nodes2 = { n3,n4 };
	std::vector<size_t> indexes2 = { 1,2 };

	Geometry geometry2(fig, fig_order, std::move(nodes2), std::move(indexes2));

	constexpr size_t axis_tag = 1;
	EXPECT_TRUE(geometry.is_periodic_pair(geometry2, axis_tag));
}

GTEST_TEST(Geometry, faces_geometry_1) {
	const Figure fig = Figure::triangle;
	const size_t fig_order = 1;

	const EuclideanVector n1 = { 1,1 };
	const EuclideanVector n2 = { 2,1 };
	const EuclideanVector n3 = { 4,2 };
	std::vector<EuclideanVector<2>> nodes = { n1,n2,n3 };
	std::vector<size_t> indexes = { 1,2,3,4 };

	Geometry geometry(fig, fig_order, std::move(nodes), std::move(indexes));
	const auto result = geometry.face_geometries();

	const Figure f_fig = Figure::line;

	const EuclideanVector f1_n1 = { 1,1 };
	const EuclideanVector f1_n2 = { 2,1 };
	std::vector<EuclideanVector<2>> f1_nodes = { f1_n1,f1_n2 };
	std::vector<size_t> f1_indexes = { 1,2 };
	Geometry f1_geometry(f_fig, fig_order, std::move(f1_nodes), std::move(f1_indexes));

	const EuclideanVector f2_n1 = { 2,1 };
	const EuclideanVector f2_n2 = { 4,2 };
	std::vector<EuclideanVector<2>> f2_nodes = { f2_n1,f2_n2 };
	std::vector<size_t> f2_indexes = { 1,2 };
	Geometry f2_geometry(f_fig, fig_order, std::move(f2_nodes), std::move(f2_indexes));

	const EuclideanVector f3_n1 = { 4,2 };
	const EuclideanVector f3_n2 = { 1,1 };
	std::vector<EuclideanVector<2>> f3_nodes = { f3_n1,f3_n2 };
	std::vector<size_t> f3_indexes = { 1,2 };
	Geometry f3_geometry(f_fig, fig_order, std::move(f3_nodes), std::move(f3_indexes));


	const std::vector<Geometry<2>> ref = { f1_geometry,f2_geometry,f3_geometry };
	EXPECT_EQ(result, ref);
}