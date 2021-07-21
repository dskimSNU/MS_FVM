#pragma once
#include "Grid_Builder.h"
#include "Governing_Equation.h"
#include "Inital_Condition.h"

class PostAI
{
//private:
public:
	inline static std::vector<std::vector<size_t>> vertex_share_cell_indexes_set_;

	inline static std::vector<std::string> node_number_string_set_;
	inline static std::vector<std::string> edge_number_string_set_;
	inline static std::vector<std::string> connectivity_string_set_;
	inline static std::vector<std::string> cell_coords_string_set_;

public:
	template <size_t space_dimension>
	static void intialize(const Grid<space_dimension>& grid);

//private: //for test
	template <size_t space_dimension>
	static auto calculate_face_share_cell_indexes_set(const Grid<space_dimension>& grid);

	template <size_t space_dimension>
	static auto calculate_vertex_nodes_coordinate_string_set(const Grid<space_dimension>& grid);
};


//template definition part
template <size_t space_dimension>
void PostAI::intialize(const Grid<space_dimension>& grid) {
	const auto& vnode_index_to_share_cell_indexes = grid.connectivity.vnode_index_to_share_cell_indexes;
	const auto& cell_elements = grid.elements.cell_elements;
	const auto num_cell = cell_elements.size();


	vertex_share_cell_indexes_set_.reserve(num_cell);
	node_number_string_set_.reserve(num_cell);
	edge_number_string_set_.reserve(num_cell);
	connectivity_string_set_.reserve(num_cell);
	cell_coords_string_set_.reserve(num_cell);


	const auto face_share_cell_indexes_set = calculate_face_share_cell_indexes_set(grid);
	const auto vnodes_coordinate_string_set = calculate_vertex_nodes_coordinate_string_set(grid);

	for (size_t i = 0; i < num_cell; ++i) {
		const auto& cell_element = cell_elements[i];
		const auto& cell_geometry = cell_element.geometry_;

		//vertex share cell indexes
		std::set<size_t> vertex_share_cell_indexes_temp;

		const auto vnode_indexes = cell_element.vertex_node_indexes();
		for (const auto& vnode_index : vnode_indexes) {
			const auto& vnode_share_cell_indexes = vnode_index_to_share_cell_indexes.at(vnode_index);
			vertex_share_cell_indexes_temp.insert(vnode_share_cell_indexes.begin(), vnode_share_cell_indexes.end());
		}

		// node number string
		const auto num_node = vertex_share_cell_indexes.size();
		node_number_string_set_.push_back("@nodeNumber\n" + std::to_string(num_node) + "\n");

		//chunk edge connectivities //quad3에서는 제대로 작동하지 않는 algorithm
		std::set<std::set<size_t>> chunk_edge_connectivities;

		for (const auto chunk_cell_index : vertex_share_cell_indexes) {
			const auto& face_share_cell_indexes = face_share_cell_indexes_set.at(chunk_cell_index);

			std::vector<size_t> face_share_cell_indexes_in_chunk;
			std::set_intersection(vertex_share_cell_indexes.begin(), vertex_share_cell_indexes.end(), face_share_cell_indexes.begin(), face_share_cell_indexes.end(), std::back_inserter(face_share_cell_indexes_in_chunk));

			for (const auto face_share_cell_index_in_chunk : face_share_cell_indexes_in_chunk)
				chunk_edge_connectivities.insert({ chunk_cell_index, face_share_cell_index_in_chunk });
		}

		//edge number string
		const auto num_edge = chunk_edge_connectivities.size();
		edge_number_string_set_.push_back("@edgeNumber\n" + std::to_string(num_edge) + "\n");

		//connectivity string
		std::string node_connectivity_string = "@connectivity\n";
		for (const auto& chunk_edge_connectivity : chunk_edge_connectivities) {
			for (const auto& node_index : chunk_edge_connectivity)
				node_connectivity_string += std::to_string(node_index) + "\t";
			node_connectivity_string += "\n";
		}
		connectivity_string_set_.push_back(std::move(node_connectivity_string));

		//cell coords string
		std::string cell_coords_string = "@cellCoords\n";
		for (const auto vertex_share_cell_index : vertex_share_cell_indexes) 
			cell_coords_string += vnodes_coordinate_string_set[vertex_share_cell_index];

		cell_coords_string_set_.push_back(std::move(cell_coords_string));

		//last
		vertex_share_cell_indexes_temp.erase(i);

		std::vector<size_t> vertex_share_cell_indexes;
		vertex_share_cell_indexes.push_back(i);
		vertex_share_cell_indexes.insert(vertex_share_cell_indexes.end(), vertex_share_cell_indexes_temp.begin(), vertex_share_cell_indexes_temp.end());

		vertex_share_cell_indexes_set_.push_back(vertex_share_cell_indexes);

	}
}


template <size_t space_dimension>
auto PostAI::calculate_face_share_cell_indexes_set(const Grid<space_dimension>& grid) {
	const auto& vnode_index_to_share_cell_indexes = grid.connectivity.vnode_index_to_share_cell_indexes;
	const auto& cell_elements = grid.elements.cell_elements;
	const auto num_cell = cell_elements.size();

	//face share cell indexes set
	std::vector<std::set<size_t>> face_share_cell_indexes_set;
	face_share_cell_indexes_set.reserve(num_cell);

	for (size_t i = 0; i < num_cell; ++i) {
		const auto& element = cell_elements[i];
		const auto& geometry = cell_elements[i].geometry_;

		const auto face_vnode_indexes_set = element.face_vertex_node_indexes_set();
		const auto num_face = face_vnode_indexes_set.size();

		std::set<size_t> face_share_cell_indexes;

		for (const auto& face_vnode_indexes : face_vnode_indexes_set) {
			std::vector<size_t> this_face_share_cell_indexes;

			const auto num_face_vnode = face_vnode_indexes.size();

			const auto& set_0 = vnode_index_to_share_cell_indexes.at(face_vnode_indexes[0]);
			const auto& set_1 = vnode_index_to_share_cell_indexes.at(face_vnode_indexes[1]);
			std::set_intersection(set_0.begin(), set_0.end(), set_1.begin(), set_1.end(), std::back_inserter(this_face_share_cell_indexes));

			if (2 < num_face_vnode) {
				std::vector<size_t> buffer;
				for (size_t i = 2; i < num_face_vnode; ++i) {
					const auto& set_i = vnode_index_to_share_cell_indexes.at(face_vnode_indexes[i]);

					buffer.clear();
					std::set_intersection(this_face_share_cell_indexes.begin(), this_face_share_cell_indexes.end(), set_i.begin(), set_i.end(), std::back_inserter(buffer));
					std::swap(this_face_share_cell_indexes, buffer);
				}
			}

			const auto my_index_pos_iter = std::find(this_face_share_cell_indexes.begin(), this_face_share_cell_indexes.end(), i);
			dynamic_require(my_index_pos_iter != this_face_share_cell_indexes.end(), "my index should be included in this face share cell indexes");

			this_face_share_cell_indexes.erase(my_index_pos_iter);
			dynamic_require(this_face_share_cell_indexes.size() == 1, "face share cell should be unique");

			face_share_cell_indexes.insert(this_face_share_cell_indexes.front());
		}

		face_share_cell_indexes_set.push_back(std::move(face_share_cell_indexes));
	}

	return face_share_cell_indexes_set;
}


template <size_t space_dimension>
auto PostAI::calculate_vertex_nodes_coordinate_string_set(const Grid<space_dimension>& grid) {
	const auto& cell_elements = grid.elements.cell_elements;
	const auto num_cell = cell_elements.size();

	std::vector<std::string> vnodes_coordinate_string_set;
	vnodes_coordinate_string_set.reserve(num_cell);

	for (size_t i = 0; i < num_cell; ++i) {
		const auto& element = cell_elements[i];
		const auto& geometry = element.geometry_;

		const auto vnodes = geometry.vertex_nodes();
		const auto num_vnode = vnodes.size();
		std::string vnodes_coordinate_string = std::to_string(num_vnode) + "\n";

		for (const auto& vnode : vnodes) {
			std::ostringstream stream;
			stream << std::setprecision(16) << std::showpoint;

			for (size_t i = 0; i < space_dimension; ++i)
				stream << vnode[i] << "\t";

			vnodes_coordinate_string += stream.str() + "\n";
		}

		vnodes_coordinate_string_set.push_back(std::move(vnodes_coordinate_string));
	}

	return vnodes_coordinate_string_set;
}
