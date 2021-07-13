#pragma once
#include "Matrix.h"
#include <set>

enum class Figure
{
	point,
	line,
	triangle, quadrilateral,
	tetrahedral, hexahedral, prism, pyramid,
	not_in_list
};


template <size_t space_dimension>
class Geometry
{
	class ReferenceGeometry;	//罹晦縑 氈戲棲梱 Test陛 ��虛! 疏頂朝啪 釭擊梱 ?

	using Space_Vector = EuclideanVector<space_dimension>;
		
public:
	Geometry(const ReferenceGeometry reference_geometry, std::vector<Space_Vector>&& consisting_nodes, std::vector<size_t>&& consisting_node_indexes)
		: reference_geometry_(reference_geometry), nodes_(std::move(consisting_nodes)), node_indexes_(std::move(consisting_node_indexes)) {};

	Geometry(const Figure figure, const size_t figure_order, std::vector<Space_Vector>&& consisting_nodes, std::vector<size_t>&& consisting_node_indexes)
		: reference_geometry_(figure, figure_order), nodes_(std::move(consisting_nodes)), node_indexes_(std::move(consisting_node_indexes)) {};


	bool operator==(const Geometry& other) const;
	bool operator<(const Geometry& other) const;

	Space_Vector center_node(void) const;
	Space_Vector normal_vector(const Space_Vector& owner_cell_center) const;
	double volume(void) const;
	std::array<double, space_dimension> coordinate_projected_volume(void) const;
	std::vector<Geometry> face_geometries(void) const;
	std::vector<Space_Vector> vertex_nodes(void) const;
	std::vector<size_t> vertex_node_indexes(void) const;
	bool is_periodic_pair(const Geometry& other, const size_t axis_tag) const;
	std::vector<std::vector<size_t>> local_connectivities(void) const;

	//private: for test
	std::vector<std::vector<Space_Vector>> calculate_faces_nodes(void) const;
	bool is_perioidc_pair_node(const Space_Vector& node, const size_t axis_tag) const;
	bool is_axis_translation(const Space_Vector& node1, const Space_Vector& node2, const size_t axis_tag) const;

private:
	ReferenceGeometry reference_geometry_;
	std::vector<Space_Vector> nodes_;
	std::vector<size_t> node_indexes_;
};


//class Reference Geometry
template <size_t space_dimension>
class Geometry<space_dimension>::ReferenceGeometry
{
public:
	ReferenceGeometry(const Figure figure, const size_t figure_order) : figure_(figure), figure_order_(figure_order) {};

	bool operator==(const ReferenceGeometry& other) const;
	bool operator!=(const ReferenceGeometry& other) const;

	std::vector<size_t> vertex_node_index_orders(void) const;
	std::vector<std::vector<size_t>> calculate_faces_node_index_orders(void) const;
	std::vector<ReferenceGeometry> calculate_faces_reference_geometry(void) const;
	Space_Vector calculate_normal(const std::vector<Space_Vector>& nodes) const;
	double calculate_volume(const std::vector<Space_Vector>& nodes) const;
	std::vector<std::vector<size_t>> local_connectivities(void) const;

private:
	Figure figure_;
	size_t figure_order_;
};


//Reference Geometry template definition part
template <size_t space_dimension>
bool Geometry<space_dimension>::ReferenceGeometry::operator==(const ReferenceGeometry& other) const {
	return this->figure_ == other.figure_ && this->figure_order_ == other.figure_order_;
}

template <size_t space_dimension>
bool Geometry<space_dimension>::ReferenceGeometry::operator != (const ReferenceGeometry & other) const {
	return !((*this) == other);
}

template <size_t space_dimension>
std::vector<size_t> Geometry<space_dimension>::ReferenceGeometry::vertex_node_index_orders(void) const {
	switch (this->figure_) {
	case Figure::line: {
		// 0 式式式式 1
		return { 0,1 };
	}
	case Figure::triangle: {
		//  2
		//  弛 \
		//	弛  \
		//  0式式式1
		return { 0,1,2 };
	}
	case Figure::quadrilateral: {
		//  3式式式式式2
		//  弛     弛
		//  0式式式式式1
		return { 0,1,2,3 };
	}
	default:
		throw std::runtime_error("wrong element figure");
		return { 0 };
	}
}

template <size_t space_dimension>
std::vector<std::vector<size_t>> Geometry<space_dimension>::ReferenceGeometry::calculate_faces_node_index_orders(void) const {
	switch (this->figure_) {
	case Figure::line: {
		// 0 式式式式 1
		const std::vector<size_t> face0_node_index = { 0 };
		const std::vector<size_t> face1_node_index = { 1 };
		return { face0_node_index,face1_node_index };
	}
	case Figure::triangle: {
		//      2
		//  2  / \  1
		//	  /   \
		//   0式式式式式1
		//      0
		const std::vector<size_t> face0_node_index = { 0,1 };
		const std::vector<size_t> face1_node_index = { 1,2 };
		const std::vector<size_t> face2_node_index = { 2,0 };
		return { face0_node_index,face1_node_index, face2_node_index };

		//if (element_order > 1)
		//{
		//	const size_t num_additional_point = element_order - 1;
		//	
		//	size_t index = num_face;
		//	for (size_t iface = 0; iface < num_face; ++iface)
		//		for (size_t ipoint = 0; ipoint < num_additional_point; ++ipoint)
		//			face_node_index_order[iface].push_back(index++);
		//}
	}
	case Figure::quadrilateral: {
		//      2
		//   3式式式式式2
		//3  弛     弛   1
		//   0式式式式式1
		//      0
		std::vector<size_t> face0_node_index = { 0,1 };
		std::vector<size_t> face1_node_index = { 1,2 };
		std::vector<size_t> face2_node_index = { 2,3 };
		std::vector<size_t> face3_node_index = { 3,0 };
		return { face0_node_index,face1_node_index, face2_node_index,face3_node_index };

		//if (element_order > 1)
		//{
		//	const size_t num_additional_point = element_order - 1;

		//	size_t index = num_face;
		//	for (size_t iface = 0; iface < num_face; ++iface)
		//		for (size_t ipoint = 0; ipoint<num_additional_point; ++ipoint)
		//			face_node_index_order[iface].push_back(index++);
		//}
	}
	default:
		throw std::runtime_error("wrong element figure");
		return std::vector<std::vector<size_t>>();
	}
}

template <size_t space_dimension>
std::vector<typename Geometry<space_dimension>::ReferenceGeometry> Geometry<space_dimension>::ReferenceGeometry::calculate_faces_reference_geometry(void) const {
	switch (this->figure_) {
	case Figure::line: {
		// 0 式式式式 1
		const ReferenceGeometry face0_reference_geometry = { Figure::point,this->figure_order_ };
		const ReferenceGeometry face1_reference_geometry = { Figure::point,this->figure_order_ };
		return { face0_reference_geometry,face1_reference_geometry };
	}
	case Figure::triangle: {
		//      2
		//  2  / \  1
		//	  /   \
		//   0式式式式式1
		//      0
		const ReferenceGeometry face0_refrence_geometry = { Figure::line,this->figure_order_ };
		const ReferenceGeometry face1_refrence_geometry = { Figure::line,this->figure_order_ };
		const ReferenceGeometry face2_refrence_geometry = { Figure::line,this->figure_order_ };
		return { face0_refrence_geometry,face1_refrence_geometry, face2_refrence_geometry };
	}
	case Figure::quadrilateral: {
		//      2
		//   3式式式式式2
		//3  弛     弛   1
		//   0式式式式式1
		//      0
		const ReferenceGeometry face0_refrence_geometry = { Figure::line,this->figure_order_ };
		const ReferenceGeometry face1_refrence_geometry = { Figure::line,this->figure_order_ };
		const ReferenceGeometry face2_refrence_geometry = { Figure::line,this->figure_order_ };
		const ReferenceGeometry face3_refrence_geometry = { Figure::line,this->figure_order_ };
		return { face0_refrence_geometry,face1_refrence_geometry, face2_refrence_geometry,face3_refrence_geometry };
	}
	default:
		throw std::runtime_error("not supported figure");
		return std::vector<ReferenceGeometry>();
	}
}

template <size_t space_dimension>
Geometry<space_dimension>::Space_Vector Geometry<space_dimension>::ReferenceGeometry::calculate_normal(const std::vector<Space_Vector>& nodes) const {
	switch (this->figure_) {
	case Figure::line: {
		// 0 式式式> 1
		const static Matrix<2, 2> rotation_matrix = { 0, -1, 1, 0 };
		const auto a = nodes[1] - nodes[0];
		return (rotation_matrix * a).be_normalize();
	}
	default:
		throw std::runtime_error("not supported element figure");
		return {};
	}
}

template <size_t space_dimension>
double Geometry<space_dimension>::ReferenceGeometry::calculate_volume(const std::vector<Space_Vector>& nodes) const {
	switch (this->figure_) {
	case Figure::line: {
		return (nodes[1] - nodes[0]).norm();
	}
	case Figure::triangle: {
		//  2
		//  ^ \ 
		//b 弛  \
		//  0式式>1
		//    a
		const auto a = nodes[1] - nodes[0];
		const auto b = nodes[2] - nodes[0];
		return 0.5 * std::sqrt(a.inner_product(a) * b.inner_product(b) - a.inner_product(b) * a.inner_product(b));

		//Heron's formula - severe round off error
		//const auto a = (*nodes[1] - *nodes[0]).norm();
		//const auto b = (*nodes[2] - *nodes[0]).norm();
		//const auto c = (*nodes[2] - *nodes[1]).norm();
		//const auto s = (a + b + c) * 0.5;
		//return std::sqrt(s * (s - a) * (s - b) * (s - c));
	}
	case Figure::quadrilateral: {
		//     c
		//  3式式式式>2
		//  弛     ^
		//d v     弛 b
		//  0<式式式式1
		//	   a

		const auto a = nodes[0] - nodes[1];
		const auto b = nodes[2] - nodes[1];
		const auto c = nodes[2] - nodes[3];
		const auto d = nodes[0] - nodes[3];

		const auto triangle_ab = 0.5 * std::sqrt(a.inner_product(a) * b.inner_product(b) - a.inner_product(b) * a.inner_product(b));
		const auto triangle_cd = 0.5 * std::sqrt(c.inner_product(c) * d.inner_product(d) - c.inner_product(d) * c.inner_product(d));

		return triangle_ab + triangle_cd;
	}
	default:
		throw std::runtime_error("wrong element figure");
		return { 0 };
	}
}

template <size_t space_dimension>
std::vector<std::vector<size_t>> Geometry<space_dimension>::ReferenceGeometry::local_connectivities(void) const {
	switch (this->figure_) {
	case Figure::triangle: {
		//  2
		//  弛 \ 
		//  弛  \
		//  0式式式1
		return { { 0,1,2 } };
	}
	case Figure::quadrilateral: {
		//  3式式式式式2
		//  弛     弛
		//  弛     弛 
		//  0式式式式式1

		return { {0,1,2},{0,2,3} };
	}
	default:
		throw std::runtime_error("wrong element figure");
		return {};
	}
}

//Geomety template definition part
template <size_t space_dimension>
std::array<double, space_dimension> Geometry<space_dimension>::coordinate_projected_volume(void) const {
	if constexpr (space_dimension == 2) {
		double x_projected_volume = 0.0;
		double y_projected_volume = 0.0;

		const auto faces_nodes = this->calculate_faces_nodes();
		for (const auto& face_nodes : faces_nodes) {
			const auto& start_node = face_nodes[0];
			const auto& end_node = face_nodes[1];
			const auto node_to_node = end_node - start_node;

			x_projected_volume += std::abs(node_to_node[0]);
			y_projected_volume += std::abs(node_to_node[1]);
		}

		return { 0.5 * x_projected_volume, 0.5 * y_projected_volume };
	}
	else {
		throw std::runtime_error("not supproted dimension");
		return {};
	}
}

template<size_t space_dimension>
std::vector<Geometry<space_dimension>> Geometry<space_dimension>::face_geometries(void) const {
	const auto faces_reference_geometry = this->reference_geometry_.calculate_faces_reference_geometry();

	const auto faces_node_index_orders = this->reference_geometry_.calculate_faces_node_index_orders();
	const auto num_face = faces_node_index_orders.size();

	std::vector<Geometry> faces_geometry;
	faces_geometry.reserve(num_face);

	for (size_t i = 0; i < num_face; ++i) {
		const auto& reference_geometry = faces_reference_geometry[i];

		const auto& node_index_orders = faces_node_index_orders[i];
		const auto num_node = node_index_orders.size();

		std::vector<Space_Vector> face_nodes(num_node);
		std::vector<size_t> face_node_indexes(num_node);
		for (size_t j = 0; j < num_node; ++j) {
			face_nodes[j] = this->nodes_[node_index_orders[j]];
			face_node_indexes[j] = this->node_indexes_[node_index_orders[j]];
		}

		faces_geometry.push_back({ faces_reference_geometry[i],std::move(face_nodes), std::move(face_node_indexes) });
	}

	return faces_geometry;
}

template<size_t space_dimension>
std::vector<EuclideanVector<space_dimension>> Geometry<space_dimension>::vertex_nodes(void) const {
	const auto vertex_node_index_orders = this->reference_geometry_.vertex_node_index_orders();
	const auto num_vertex_node = vertex_node_index_orders.size();

	std::vector<Space_Vector> vertex_nodes(num_vertex_node);
	for (size_t i = 0; i < num_vertex_node; ++i)
		vertex_nodes[i] = this->nodes_[vertex_node_index_orders[i]];

	return vertex_nodes;
}

template<size_t space_dimension>
std::vector<size_t> Geometry<space_dimension>::vertex_node_indexes(void) const {
	const auto vertex_node_index_orders = this->reference_geometry_.vertex_node_index_orders();
	const auto num_vertex_node = vertex_node_index_orders.size();

	std::vector<size_t> vertex_node_indexes(num_vertex_node);
	for (size_t i = 0; i < num_vertex_node; ++i)
		vertex_node_indexes[i] = this->node_indexes_[vertex_node_index_orders[i]];

	return vertex_node_indexes;
}

template<size_t space_dimension>
std::vector<std::vector<typename Geometry<space_dimension>::Space_Vector>> Geometry<space_dimension>::calculate_faces_nodes(void) const {
	const auto faces_node_index_orders = this->reference_geometry_.calculate_faces_node_index_orders();
	const auto num_face = faces_node_index_orders.size();

	std::vector<std::vector<Space_Vector>> faces_nodes(num_face);
	for (size_t i = 0; i < num_face; ++i) {
		const auto& face_node_index_orders = faces_node_index_orders[i];
		const auto num_node = face_node_index_orders.size();

		auto& face_nodes = faces_nodes[i];
		face_nodes.resize(num_node);

		for (size_t j = 0; j < face_node_index_orders.size(); ++j)
			face_nodes[j] = this->nodes_[face_node_index_orders[j]];
	}

	return faces_nodes;
}

template<size_t space_dimension>
bool Geometry<space_dimension>::is_periodic_pair(const Geometry& other, const size_t axis_tag) const {
	if (this->reference_geometry_ != other.reference_geometry_)
		return false;

	for (const auto& node : other.nodes_) {
		if (this->is_perioidc_pair_node(node, axis_tag))
			continue;
		else
			return false;
	}
	return true;
}

template <size_t space_dimension>
std::vector<std::vector<size_t>> Geometry<space_dimension>::local_connectivities(void) const {
	return this->reference_geometry_.local_connectivities();
}

template<size_t space_dimension>
bool Geometry<space_dimension>::is_perioidc_pair_node(const Space_Vector& node, const size_t axis_tag) const {
	for (const auto& my_node : this->nodes_) {
		if (this->is_axis_translation(my_node, node, axis_tag))
			return true;
	}
	return false;
}

template<size_t space_dimension>
bool Geometry<space_dimension>::is_axis_translation(const Space_Vector& node1, const Space_Vector& node2, const size_t axis_tag) const {
	const auto line_vector = node2 - node1;
	for (size_t i = 0; i < space_dimension; ++i) {
		if (i == axis_tag)
			continue;

		if (std::abs(line_vector[i]) > 1.0E-10)
			return false;
	}
	return true;
}

template<size_t space_dimension>
bool Geometry<space_dimension>::operator==(const Geometry& other) const {
	return this->reference_geometry_ == other.reference_geometry_ && this->nodes_ == other.nodes_;
}

template<size_t space_dimension>
bool Geometry<space_dimension>::operator<(const Geometry& other) const {
	std::set<size_t> s1(this->node_indexes_.begin(), this->node_indexes_.end());
	std::set<size_t> s2(other.node_indexes_.begin(), other.node_indexes_.end());
	return s1 > s2;
}

template<size_t space_dimension>
Geometry<space_dimension>::Space_Vector Geometry<space_dimension>::center_node(void) const {
	Space_Vector center;
	for (const auto& node : this->nodes_)
		center += node;

	const auto num_node = this->nodes_.size();
	return center * (1.0 / num_node);
}

template<size_t space_dimension>
Geometry<space_dimension>::Space_Vector Geometry<space_dimension>::normal_vector(const Space_Vector& owner_cell_center) const {
	const auto normal = this->reference_geometry_.calculate_normal(this->nodes_);
	const auto vector_pointing_outward = this->center_node() - owner_cell_center;

	if (normal.inner_product(vector_pointing_outward) > 0)
		return normal;
	else
		return -1 * normal;

}

template<size_t space_dimension>
double Geometry<space_dimension>::volume(void) const {
	return this->reference_geometry_.calculate_volume(this->nodes_);
}