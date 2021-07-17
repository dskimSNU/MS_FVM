#pragma once
#include "Matrix.h"
#include <set>
#include <unordered_set>

enum class Figure
{
	point,
	line,
	triangle, quadrilateral,
	tetrahedral, hexahedral, prism, pyramid,
	not_in_list
};


enum class ElementType
{
	cell, inner_face,
	slip_wall_2D,
	supersonic_outlet_2D,
	x_periodic, y_periodic,
	not_in_list
};


class ReferenceGeometry
{
public:
	ReferenceGeometry(const Figure figure, const size_t figure_order) : figure_(figure), figure_order_(figure_order) {};

	bool operator==(const ReferenceGeometry& other) const;
	bool operator!=(const ReferenceGeometry& other) const;

	std::vector<size_t> vertex_node_index_orders(void) const;
	std::vector<std::vector<size_t>> faces_node_index_orders(void) const;
	std::vector<ReferenceGeometry> faces_reference_geometry(void) const;
	std::vector<std::vector<size_t>> local_connectivities(void) const;

	template <typename Space_Vector>
	Space_Vector calculate_normal(const std::vector<Space_Vector>& nodes) const;
	template <typename Space_Vector>
	double calculate_volume(const std::vector<Space_Vector>& nodes) const;


private:
	Figure figure_;
	size_t figure_order_;
};


template <size_t space_dimension>
class Geometry
{
private:
	using Space_Vector = EuclideanVector<space_dimension>;

public:
	ReferenceGeometry reference_geometry_;

private:	
	std::vector<Space_Vector> nodes_;
		
public:
	Geometry(const ReferenceGeometry reference_geometry, std::vector<Space_Vector>&& consisting_nodes)
		: reference_geometry_(reference_geometry), nodes_(std::move(consisting_nodes)) {};

	Space_Vector center_node(void) const;
	Space_Vector normal_vector(const Space_Vector& owner_cell_center) const;
	double volume(void) const;
	std::array<double, space_dimension> coordinate_projected_volume(void) const;
	std::vector<Geometry> faces_geometry(void) const;
	std::vector<Space_Vector> vertex_nodes(void) const;
	bool is_axis_parallel(const Geometry& other, const size_t axis_tag) const;

	//private: for test
	std::vector<std::vector<Space_Vector>> calculate_faces_nodes(void) const;
	bool is_axis_parallel_node(const Space_Vector& node, const size_t axis_tag) const;
};


template <size_t space_dimension>
class Element
{
public:
	Geometry<space_dimension> geometry_;

private:
	ElementType element_type_;
	std::vector<size_t> node_indexes_;

public:
	Element(const ElementType element_type, Geometry<space_dimension>&& geometry, std::vector<size_t>&& node_indexes)
		: element_type_(element_type), geometry_(std::move(geometry)), node_indexes_(std::move(node_indexes)) {};

	std::vector<size_t> vertex_node_indexes(void) const;
	std::vector<Element> make_inner_face_elements(void) const;
	bool is_periodic_pair(const Element& other) const;
	std::vector<std::pair<size_t, size_t>> find_periodic_vnode_index_pairs(const Element& other) const;

//private:
	std::vector<std::vector<size_t>> faces_node_indexes(void) const;
	bool is_periodic_boundary(void) const;

};


//template definition part

template <typename Space_Vector>
Space_Vector ReferenceGeometry::calculate_normal(const std::vector<Space_Vector>& nodes) const {
	switch (this->figure_) {
	case Figure::line: {
		// 0 式式式> 1
		//    a
		dynamic_require(nodes.size() == 2, "line should have two vertex");
		const static Matrix<2, 2> rotation_matrix = { 0, -1, 1, 0 };
		const auto a = nodes[1] - nodes[0];
		return (rotation_matrix * a).be_normalize();
	}
	default:
		throw std::runtime_error("not supported element figure");
		return {};
	}
}

template <typename Space_Vector>
double ReferenceGeometry::calculate_volume(const std::vector<Space_Vector>& nodes) const {
	switch (this->figure_) {
	case Figure::line: {
		// 0 式式式> 1
		dynamic_require(nodes.size() == 2, "line should have two vertex");
		return (nodes[1] - nodes[0]).norm();
	}
	case Figure::triangle: {
		//  2
		//  ^ \ 
		//b 弛  \
		//  0式式>1
		//    a
		dynamic_require(nodes.size() == 3, "triangle should have three vertex");

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
		dynamic_require(nodes.size() == 4, "quadrilateral should have four vertex");

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
std::vector<Geometry<space_dimension>> Geometry<space_dimension>::faces_geometry(void) const {
	const auto faces_reference_geometry = this->reference_geometry_.faces_reference_geometry();
	const auto faces_node_index_orders = this->reference_geometry_.faces_node_index_orders();
	const auto num_face = faces_node_index_orders.size();

	std::vector<Geometry> faces_geometry;
	faces_geometry.reserve(num_face);

	for (size_t i = 0; i < num_face; ++i) {
		const auto& reference_geometry = faces_reference_geometry[i];

		const auto& node_index_orders = faces_node_index_orders[i];
		const auto num_node = node_index_orders.size();

		std::vector<Space_Vector> face_nodes(num_node);
		for (size_t j = 0; j < num_node; ++j) 
			face_nodes[j] = this->nodes_[node_index_orders[j]];

		faces_geometry.push_back({ reference_geometry,std::move(face_nodes) });
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
std::vector<std::vector<typename Geometry<space_dimension>::Space_Vector>> Geometry<space_dimension>::calculate_faces_nodes(void) const {
	const auto faces_node_index_orders = this->reference_geometry_.faces_node_index_orders();
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
bool Geometry<space_dimension>::is_axis_parallel(const Geometry& other, const size_t axis_tag) const {
	if (this->reference_geometry_ != other.reference_geometry_)
		return false;

	if (this->nodes_.size() != other.nodes_.size())
		return false;

	for (const auto& node : other.nodes_) {
		if (this->is_axis_parallel_node(node, axis_tag))
			continue;
		else
			return false;
	}
	return true;
}

template<size_t space_dimension>
bool Geometry<space_dimension>::is_axis_parallel_node(const Space_Vector& node, const size_t axis_tag) const {
	for (const auto& my_node : this->nodes_) {
		if (my_node.is_axis_translation(node,axis_tag))
			return true;
	}
	return false;
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

template<size_t space_dimension>
std::vector<size_t> Element<space_dimension>::vertex_node_indexes(void) const {
	const auto vertex_node_index_orders = this->geometry_.reference_geometry_.vertex_node_index_orders();
	const auto num_vertex_node = vertex_node_index_orders.size();

	std::vector<size_t> vertex_node_indexes(num_vertex_node);
	for (size_t i = 0; i < num_vertex_node; ++i)
		vertex_node_indexes[i] = this->node_indexes_[vertex_node_index_orders[i]];

	return vertex_node_indexes;
}

template<size_t space_dimension>
std::vector<Element<space_dimension>> Element<space_dimension>::make_inner_face_elements(void) const {
	dynamic_require(this->element_type_ == ElementType::cell, "make inner face elements should be called from cell element");

	auto faces_geometry = this->geometry_.faces_geometry();
	auto faces_node_indexes = this->faces_node_indexes();

	const auto num_face = faces_geometry.size();
	std::vector<Element<space_dimension>> inner_face_elements;
	inner_face_elements.reserve(num_face);

	for (size_t i = 0; i < num_face; ++i) 
		inner_face_elements.push_back({ ElementType::inner_face, std::move(faces_geometry[i]), std::move(faces_node_indexes[i]) });

	return inner_face_elements;
}


template<size_t space_dimension>
std::vector<std::vector<size_t>> Element<space_dimension>::faces_node_indexes(void) const {
	const auto faces_node_index_orders = this->geometry_.reference_geometry_.faces_node_index_orders();
	const auto num_face = faces_node_index_orders.size();

	std::vector<std::vector<size_t>> faces_node_indexes(num_face);
	for (size_t i = 0; i < num_face; ++i) {
		const auto& face_node_index_orders = faces_node_index_orders[i];
		const auto num_node = face_node_index_orders.size();

		auto& face_node_indexes = faces_node_indexes[i];
		face_node_indexes.resize(num_node);

		for (size_t j = 0; j < num_node; ++j)
			face_node_indexes[j] = this->node_indexes_[face_node_index_orders[j]];
	}

	return faces_node_indexes;
}

template<size_t space_dimension>
bool Element<space_dimension>::is_periodic_pair(const Element& other) const {
	dynamic_require(this->is_periodic_boundary() && other.is_periodic_boundary(), "both elemets should be periodic boundary");

	if (this->element_type_ != other.element_type_)
		return false;

	size_t axis_tag;
	if (this->element_type_ == ElementType::x_periodic)
		axis_tag = 0;
	else
		axis_tag = 1;

	if (this->geometry_.is_axis_parallel(other.geometry_, axis_tag))
		return true;
	else
		return false;
}

template<size_t space_dimension>
std::vector<std::pair<size_t, size_t>> Element<space_dimension>::find_periodic_vnode_index_pairs(const Element& other) const {
	//both element should be periodic pair
	const auto this_vnode_indexes = this->vertex_node_indexes();
	const auto other_vnode_indexes = other.vertex_node_indexes();

	const auto this_vnodes = this->geometry_.vertex_nodes();
	const auto other_vnodes = other.geometry_.vertex_nodes();

	const auto this_num_vnode = this_vnode_indexes.size();
	const auto other_num_vnode = other_vnode_indexes.size();
	dynamic_require(this_num_vnode == other_num_vnode, "periodic pair should have same number of vertex node");
	
	size_t axis_tag = 0;
	if (this->element_type_ == ElementType::x_periodic)
		axis_tag = 0;
	else if (this->element_type_ == ElementType::y_periodic)
		axis_tag = 1;
	else
		throw std::runtime_error("wrong element type");
	dynamic_require(this->element_type_ == other.element_type_, "periodic pair should have same element type");

	std::unordered_set<size_t> matched_other_vnode_index;
	matched_other_vnode_index.reserve(other_num_vnode);

	std::vector<std::pair<size_t, size_t>> periodic_vnode_index_pairs;
	periodic_vnode_index_pairs.reserve(this_num_vnode);

	for (size_t i = 0; i < this_num_vnode; ++i) {
		const auto& this_vnode = this_vnodes[i];
		const auto this_vnode_index = this_vnode_indexes[i];
		for (size_t j = 0; j < other_num_vnode; ++j) {
			const auto& other_vnode = other_vnodes[j];
			const auto other_vnode_index = other_vnode_indexes[j];

			if (matched_other_vnode_index.find(other_vnode_index) != matched_other_vnode_index.end())
				continue;

			if (this_vnode.is_axis_translation(other_vnode, axis_tag)) {
				periodic_vnode_index_pairs.push_back(std::make_pair(this_vnode_index, other_vnode_index));
				matched_other_vnode_index.emplace(other_vnode_index);
			}
		}
	}
	
	dynamic_require(periodic_vnode_index_pairs.size() == this_num_vnode, "every vnode should have pair");
	return periodic_vnode_index_pairs;
}


template<size_t space_dimension>
bool Element<space_dimension>::is_periodic_boundary(void) const {
	switch (this->element_type_)	{
		case ElementType::x_periodic:
		case ElementType::y_periodic:	return true;
		default:						return false;
	}
}