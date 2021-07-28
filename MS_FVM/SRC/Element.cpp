#include "../INC/Element.h"

bool ReferenceGeometry::operator==(const ReferenceGeometry& other) const {
	return this->figure_ == other.figure_ && this->figure_order_ == other.figure_order_;
}

bool ReferenceGeometry::operator != (const ReferenceGeometry& other) const {
	return !((*this) == other);
}

size_t ReferenceGeometry::num_vertex(void) const {
	switch (this->figure_) {
	case Figure::line:			return 2;
	case Figure::triangle:		return 3;
	case Figure::quadrilateral:	return 4;
	default:
		throw std::runtime_error("wrong element figure");
		return NULL;
	}
}

std::vector<order> ReferenceGeometry::vertex_node_index_orders(void) const {
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

std::vector<std::vector<order>> ReferenceGeometry::face_vertex_node_index_orders_set(void) const {
	switch (this->figure_) {
	case Figure::line: {
		// 0 式式式式 1
		const std::vector<order> face0_node_index = { 0 };
		const std::vector<order> face1_node_index = { 1 };
		return { face0_node_index,face1_node_index };
	}
	case Figure::triangle: {
		//      2
		//  2  / \  1
		//	  /   \
		//   0式式式式式1
		//      0
		const std::vector<order> face0_node_index = { 0,1 };
		const std::vector<order> face1_node_index = { 1,2 };
		const std::vector<order> face2_node_index = { 2,0 };
		return { face0_node_index,face1_node_index, face2_node_index };
	}
	case Figure::quadrilateral: {
		//      2
		//   3式式式式式2
		//3  弛     弛   1
		//   0式式式式式1
		//      0
		std::vector<order> face0_node_index = { 0,1 };
		std::vector<order> face1_node_index = { 1,2 };
		std::vector<order> face2_node_index = { 2,3 };
		std::vector<order> face3_node_index = { 3,0 };
		return { face0_node_index,face1_node_index, face2_node_index,face3_node_index };
	}
	default:
		throw std::runtime_error("wrong element figure");
		return std::vector<std::vector<order>>();
	}
}

std::vector<std::vector<order>> ReferenceGeometry::face_node_index_orders_set(void) const {
	switch (this->figure_) {
	case Figure::line: {
		// 0 式式式式 1
		const std::vector<order> face0_node_index = { 0 };
		const std::vector<order> face1_node_index = { 1 };
		return { face0_node_index,face1_node_index };
	}
	case Figure::triangle: {
		//      2
		//  2  / \  1
		//	  /   \
		//   0式式式式式1
		//      0
		constexpr count num_face = 3;
		std::vector<std::vector<order>> face_node_index_orders_set(num_face);
		face_node_index_orders_set[0] = { 0,1 };
		face_node_index_orders_set[1] = { 1,2 };
		face_node_index_orders_set[2] = { 2,0 };

		if (this->figure_order_ > 1) {
			const count num_additional_point = this->figure_order_ - 1;

			index idx = num_face;
			for (index iface = 0; iface < num_face; ++iface)
				for (index ipoint = 0; ipoint < num_additional_point; ++ipoint)
					face_node_index_orders_set[iface].push_back(idx++);
		}

		return face_node_index_orders_set;
	}
	case Figure::quadrilateral: {
		//      2
		//   3式式式式式2
		//3  弛     弛   1
		//   0式式式式式1
		//      0
		constexpr count num_face = 4;
		std::vector<std::vector<order>> face_node_index_orders_set(num_face);
		face_node_index_orders_set[0] = { 0,1 };
		face_node_index_orders_set[1] = { 1,2 };
		face_node_index_orders_set[2] = { 2,3 };
		face_node_index_orders_set[3] = { 3,0 };

		if (this->figure_order_ > 1) {
			const order num_additional_point = this->figure_order_ - 1;

			index idx = num_face;
			for (index iface = 0; iface < num_face; ++iface)
				for (index ipoint = 0; ipoint < num_additional_point; ++ipoint)
					face_node_index_orders_set[iface].push_back(idx++);
		}

		return face_node_index_orders_set;

	}
	default:
		throw std::runtime_error("wrong element figure");
		return std::vector<std::vector<order>>();
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

std::vector<std::vector<order>> ReferenceGeometry::local_connectivities(void) const {
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