#pragma once
#include "Grid_Reader.h"

class ReferenceGeometry
{
public:
	ReferenceGeometry(const ElementFigure figure, const size_t figure_order) : figure_(figure), figure_order_(figure_order) {};

	std::vector<size_t> calculate_vertex_node_index_orders(void) const;
	std::vector<std::vector<size_t>> calculate_face_node_index_orders(void) const;
	double calculate_volume(const std::vector<Physical_Domain_Vector>& nodes) const;
	size_t dimension(void) const;

private:
	ElementFigure figure_;
	size_t figure_order_;
};


class Geometry
{
public:
	Geometry(const ReferenceGeometry& reference_geometry, std::vector<Physical_Domain_Vector>&& consisting_nodes)
		: reference_geometry_(reference_geometry), consisting_nodes_(std::move(consisting_nodes)) {};

	double calculate_volume(void) const;
	std::array<double, s_physical_domain_dimension> calculate_coordinate_projected_volume(void) const;
	std::vector<size_t> calculate_vertex_node_indexes(const std::vector<size_t>& consisting_node_indexes) const;

//private: for test
	std::vector<std::vector<Physical_Domain_Vector>> calculate_faces_nodes(void) const;

private:
	ReferenceGeometry reference_geometry_;
	std::vector<Physical_Domain_Vector> consisting_nodes_;
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