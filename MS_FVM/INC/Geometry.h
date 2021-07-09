#pragma once
#include "Grid_Reader.h"
#include "Matrix.h"

class ReferenceGeometry
{
public:
	ReferenceGeometry(const ElementFigure figure, const size_t figure_order) : figure_(figure), figure_order_(figure_order) {};

	bool operator==(const ReferenceGeometry& other) const;
	bool operator!=(const ReferenceGeometry& other) const;

	std::vector<size_t> vertex_node_index_orders(void) const;
	std::vector<std::vector<size_t>> calculate_faces_node_index_orders(void) const;
	std::vector<ReferenceGeometry> calculate_faces_reference_geometry(void) const;
	Physical_Domain_Vector calculate_normal(const std::vector<Physical_Domain_Vector>& nodes) const;
	double calculate_volume(const std::vector<Physical_Domain_Vector>& nodes) const;
	size_t dimension(void) const;

private:
	ElementFigure figure_;
	size_t figure_order_;
};


class Geometry
{
public:
	Geometry(const ReferenceGeometry& reference_geometry, std::vector<Physical_Domain_Vector>&& consisting_nodes, std::vector<size_t>&& consisting_node_indexes)
		: reference_geometry_(reference_geometry), nodes_(std::move(consisting_nodes)), node_indexes_(std::move(consisting_node_indexes)) {};

	bool operator==(const Geometry& other) const;
	bool operator<(const Geometry& other) const;

	Physical_Domain_Vector center_node(void) const;
	Physical_Domain_Vector normal_vector(const Physical_Domain_Vector& owner_cell_center) const;
	double volume(void) const;
	std::array<double, s_physical_domain_dimension> coordinate_projected_volume(void) const;

	//std::vector<std::vector<size_t>> calculate_faces_node_indexes(const std::vector<size_t>& consisting_node_indexes) const;
	std::vector<Geometry> face_geometries(void) const;
	std::vector<size_t> vertex_node_indexes(void) const;

	bool is_periodic_pair(const Geometry& other, const ElementType type) const;

//private: for test
	std::vector<std::vector<Physical_Domain_Vector>> calculate_faces_nodes(void) const;
	bool is_perioidc_pair_node(const Physical_Domain_Vector& node, const ElementType type) const;
	bool is_axis_translation(const Physical_Domain_Vector& node1, const Physical_Domain_Vector& node2, const size_t axis_tag) const;

private:
	ReferenceGeometry reference_geometry_;
	std::vector<Physical_Domain_Vector> nodes_;
	std::vector<size_t> node_indexes_;
};

//std::vector<size_t> calculate_vertex_node_index_set(void) const {
//	const auto vertex_node_index_order_set = ms::calculate_vertex_node_index_order_set(this->figure_);
//	const auto num_vertex = vertex_node_index_order_set.size();
//
//	std::vector<size_t> vertex_node_index_set(num_vertex);
//	for (size_t i = 0; i < num_vertex; ++i)
//		vertex_node_index_set[i] = this->node_indexes[vertex_node_index_order_set[i]];
//
//	return vertex_node_index_set;
//}