#include "INC/Element.h"

bool ReferenceGeometry::operator==(const ReferenceGeometry& other) const {
	return this->figure_ == other.figure_ && this->figure_order_ == other.figure_order_;
}

bool ReferenceGeometry::operator != (const ReferenceGeometry& other) const {
	return !((*this) == other);
}

std::vector<size_t> ReferenceGeometry::vertex_node_index_orders(void) const {
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

std::vector<std::vector<size_t>> ReferenceGeometry::faces_node_index_orders(void) const {
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

std::vector<ReferenceGeometry> ReferenceGeometry::faces_reference_geometry(void) const {
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

std::vector<std::vector<size_t>> ReferenceGeometry::local_connectivities(void) const {
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