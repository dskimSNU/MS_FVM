#include "../INC/Geometry.h"

std::vector<size_t> ReferenceGeometry::calculate_vertex_node_index_orders(void) const {
	switch (this->figure_) {
	case ElementFigure::line: {
		// 0 式式式式 1
		return { 0,1 };
	}
	case ElementFigure::triangle: {
		//  2
		//  弛 \
		//	弛  \
		//  0式式式1
		return { 0,1,2 };
	}
	case ElementFigure::quadrilateral: {
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

std::vector<std::vector<size_t>> ReferenceGeometry::calculate_face_node_index_orders(void) const {
	switch (this->figure_) {
	case ElementFigure::line: {
		// 0 式式式式 1
		const std::vector<size_t> face0_node_index = { 0 };
		const std::vector<size_t> face1_node_index = { 1 };
		return { face0_node_index,face1_node_index };
	}
	case ElementFigure::triangle: {
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
	case ElementFigure::quadrilateral: {
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

double ReferenceGeometry::calculate_volume(const std::vector<Physical_Domain_Vector>& nodes) const {
	switch (this->figure_) {
	case ElementFigure::line: {
		return (nodes[1] - nodes[0]).norm();
	}
	case ElementFigure::triangle: {
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
	case ElementFigure::quadrilateral: {
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

size_t ReferenceGeometry::dimension(void) const {
	switch (this->figure_)	{
	case ElementFigure::line:			return 1;
	case ElementFigure::triangle:
	case ElementFigure::quadrilateral:	return 2;
	case ElementFigure::tetrahedral:
	case ElementFigure::hexahedral:
	case ElementFigure::prism:
	case ElementFigure::pyramid:		return 3;
	default:
		throw std::runtime_error("wrong element figure type");
		break;
	}
}

std::array<double, s_physical_domain_dimension> Geometry::calculate_coordinate_projected_volume(void) const {
	if (this->reference_geometry_.dimension() == 2) {
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
		return std::array<double, s_physical_domain_dimension>();
	}
}

std::vector<size_t> Geometry::calculate_vertex_node_indexes(const std::vector<size_t>& consisting_node_indexes) const {	
	const auto vertex_node_index_orders = this->reference_geometry_.calculate_vertex_node_index_orders();
	const auto num_vertex_node = vertex_node_index_orders.size();

	std::vector<size_t> vertex_node_indexes(num_vertex_node);
	for (size_t i = 0; i < num_vertex_node; ++i)
		vertex_node_indexes[i] = consisting_node_indexes[vertex_node_index_orders[i]];

	return vertex_node_indexes;
}

std::vector<std::vector<Physical_Domain_Vector>> Geometry::calculate_faces_nodes(void) const {
	const auto faces_node_index_orders = this->reference_geometry_.calculate_face_node_index_orders();
	const auto num_face = faces_node_index_orders.size();

	std::vector<std::vector<Physical_Domain_Vector>> faces_nodes(num_face);
	for (size_t i = 0; i < num_face; ++i) {
		const auto& face_node_index_orders = faces_node_index_orders[i];
		const auto num_node = face_node_index_orders.size();

		auto& face_nodes = faces_nodes[i];
		face_nodes.resize(num_node);

		for (size_t j = 0; j < face_node_index_orders.size(); ++j)
			face_nodes[j] = this->consisting_nodes_[face_node_index_orders[j]];
	}

	return faces_nodes;
}

double Geometry::calculate_volume(void) const {
	return this->reference_geometry_.calculate_volume(this->consisting_nodes_);
}