#pragma once
#include "Grid_Element_Builder.h"

#include <set>
#include <unordered_map>
#include <unordered_set>


template <size_t space_dimension>
struct Grid_Connectivity
{
	std::unordered_map<size_t, std::set<size_t>> vnode_index_to_share_cell_indexes; // vnode				:= node at vertex, 
																					// share_cell_indexes	:= cell indexes which shares vnode
	std::vector<size_t> boundary_oc_indexes;
	std::vector<EuclideanVector<space_dimension>> boundary_normals;

	std::vector<std::pair<size_t, size_t>> periodic_boundary_oc_nc_index_pairs;		// {owner cell index, neighbor cell index}
	std::vector<EuclideanVector<space_dimension>> periodic_boundary_normals;

	std::vector<std::pair<size_t, size_t>> inner_face_oc_nc_index_pairs;			// {owner cell index, neighbor cell index}
	std::vector<EuclideanVector<space_dimension>> inner_face_normals;
};


template <size_t space_dimension>
struct Grid
{
	Grid_Elements<space_dimension> elements;
	Grid_Connectivity<space_dimension> connectivity;
};


template<size_t space_dimension>
class Grid_Builder
{
private:
	using Space_Vector_		= EuclideanVector<space_dimension>;

public:
	template <typename Grid_File_Type>
	static Grid<space_dimension> build(const std::string& grid_file_name);

private:
	static Grid_Connectivity<space_dimension> make_grid_connectivity(const Grid_Elements<space_dimension>& grid_elements);
	static std::vector<size_t> find_cell_indexes_have_these_vnodes(const std::unordered_map<size_t, std::set<size_t>>& vertex_node_index_to_cell_container_indexes, const std::vector<size_t>& face_node_indexes);
};


//template definition part
template <size_t space_dimension>
template <typename Grid_File_Type>
Grid<space_dimension> Grid_Builder<space_dimension>::build(const std::string& grid_file_name) {
	static_require(ms::is_grid_file_type<Grid_File_Type>, "It should be grid file type");

	const auto grid_elements		= Grid_Element_Builder<Grid_File_Type, space_dimension>::build_from_grid_file(grid_file_name);
	const auto grid_connectivity	= make_grid_connectivity(grid_elements);
	return { grid_elements,grid_connectivity };
}


template <size_t space_dimension>
Grid_Connectivity<space_dimension> Grid_Builder<space_dimension>::make_grid_connectivity(const Grid_Elements<space_dimension>& grid_elements) {
	SET_TIME_POINT;

	const auto& [cell_elements, boundary_elements, periodic_boundary_element_pairs, inner_face_elements] = grid_elements;

	std::unordered_map<size_t, std::set<size_t>> vnode_index_to_share_cell_indexes;

	const auto num_cell = cell_elements.size();
	for (size_t i = 0; i < num_cell; ++i) {
		const auto vnode_indexes = cell_elements[i].vertex_node_indexes();
		for (const auto& vnode_index : vnode_indexes) {
			if (vnode_index_to_share_cell_indexes.find(vnode_index) == vnode_index_to_share_cell_indexes.end())
				vnode_index_to_share_cell_indexes.emplace(vnode_index, std::set<size_t>());

			vnode_index_to_share_cell_indexes.at(vnode_index).insert(i);
		}
	}

	//boudnary grid connectivity
	const auto num_boundary = boundary_elements.size();
	std::vector<size_t> boudnary_oc_indexes(num_boundary);
	std::vector<Space_Vector_> boundary_normals(num_boundary);

	for (size_t i = 0; i < num_boundary; ++i) {
		const auto& boundary_element = boundary_elements[i];

		const auto vnode_indexes = boundary_element.vertex_node_indexes();
		const auto oc_indexes = find_cell_indexes_have_these_vnodes(vnode_index_to_share_cell_indexes, vnode_indexes);
		dynamic_require(oc_indexes.size() == 1, "boundary should have unique owner cell");

		const auto oc_index = oc_indexes.front();
		const auto& oc_element = cell_elements.at(oc_index);
		const auto oc_center = oc_element.geometry_.center_node();

		const auto boundary_normal = boundary_element.geometry_.normal_vector(oc_center);

		boudnary_oc_indexes[i] = oc_index;
		boundary_normals[i] = boundary_normal;
	}

	//periodic boundary grid connectivity
	const auto num_pbdry_pair = periodic_boundary_element_pairs.size();
	std::vector<std::pair<size_t, size_t>> periodic_boundary_oc_nc_index_pairs(num_pbdry_pair);
	std::vector<Space_Vector_> periodic_boundary_normals(num_pbdry_pair);

	for (size_t i = 0; i < num_pbdry_pair; ++i) {
		const auto& [i_pbdry_element, j_pbdry_element] = periodic_boundary_element_pairs[i];

		const auto cell_indexes_have_i = find_cell_indexes_have_these_vnodes(vnode_index_to_share_cell_indexes, i_pbdry_element.vertex_node_indexes());
		const auto cell_indexes_have_j = find_cell_indexes_have_these_vnodes(vnode_index_to_share_cell_indexes, j_pbdry_element.vertex_node_indexes());
		dynamic_require(cell_indexes_have_i.size() == 1, "periodic boundary should have unique owner cell");
		dynamic_require(cell_indexes_have_j.size() == 1, "periodic boundary should have unique neighbor cell");

		// set i as owner side element and j as neighbor side elemnt
		const auto oc_index = cell_indexes_have_i.front();
		const auto nc_index = cell_indexes_have_j.front();

		const auto& oc_element = cell_elements[oc_index];
		const auto oc_center = oc_element.geometry_.center_node();

		const auto pbdry_normal = i_pbdry_element.geometry_.normal_vector(oc_center);

		periodic_boundary_oc_nc_index_pairs[i] = { oc_index,nc_index };
		periodic_boundary_normals[i] = pbdry_normal;
	}

	// update vnode_index_to_share_cell_indexes
	for (size_t i = 0; i < num_pbdry_pair; ++i) {
		const auto [oc_index, nc_index] = periodic_boundary_oc_nc_index_pairs[i];
		const auto& [oc_side_element, nc_side_element] = periodic_boundary_element_pairs[i];

		const auto oc_side_vnode_indexes = oc_side_element.vertex_node_indexes();
		const auto nc_side_vnode_indexes = nc_side_element.vertex_node_indexes();

		for (const auto vnode_index : oc_side_vnode_indexes) {
			auto& vnode_share_cell_indexes = vnode_index_to_share_cell_indexes.at(vnode_index);
			vnode_share_cell_indexes.insert(nc_index);
		}
		for (const auto vnode_index : nc_side_vnode_indexes) {
			auto& cell_container_indexes = vnode_index_to_share_cell_indexes.at(vnode_index);
			cell_container_indexes.insert(oc_index);
		}
	}

	for (size_t i = 0; i < num_pbdry_pair; ++i) {
		const auto& [oc_side_element, nc_side_element] = periodic_boundary_element_pairs[i];

		const auto periodic_vnode_index_pairs = oc_side_element.find_periodic_vnode_index_pairs(nc_side_element);
		for (const auto [i_vnode_index, j_vnode_index] : periodic_vnode_index_pairs) {
			auto& i_vnode_share_cell_indexes = vnode_index_to_share_cell_indexes.at(i_vnode_index);
			auto& j_vnode_share_cell_indexes = vnode_index_to_share_cell_indexes.at(j_vnode_index);

			if (i_vnode_share_cell_indexes != j_vnode_share_cell_indexes) {
				i_vnode_share_cell_indexes.insert(j_vnode_share_cell_indexes.begin(), j_vnode_share_cell_indexes.end());
				j_vnode_share_cell_indexes.insert(i_vnode_share_cell_indexes.begin(), i_vnode_share_cell_indexes.end());
			}
		}
	}


	//inner face grid connectivity
	const auto num_inner_face = inner_face_elements.size();
	std::vector<std::pair<size_t, size_t>> inner_face_oc_nc_index_pairs(num_inner_face);
	std::vector<Space_Vector_> inner_face_normals(num_inner_face);

	for (size_t i = 0; i < num_inner_face; ++i) {
		const auto& inner_face_element = inner_face_elements[i];

		const auto cell_indexes = find_cell_indexes_have_these_vnodes(vnode_index_to_share_cell_indexes, inner_face_element.vertex_node_indexes());
		dynamic_require(cell_indexes.size() == 2, "inner face should have owner cell and neighbor cell");

		//set first index as oc index
		const auto oc_index = cell_indexes[0];
		const auto nc_index = cell_indexes[1];

		const auto& oc_element = cell_elements[oc_index];
		const auto oc_center = oc_element.geometry_.center_node();

		const auto inner_face_normal = inner_face_element.geometry_.normal_vector(oc_center);

		inner_face_oc_nc_index_pairs[i] = { oc_index,nc_index };
		inner_face_normals[i] = inner_face_normal;
	}

	Log::content_ << std::left << std::setw(50) << "@ Figure out connectivity" << " ----------- " << GET_TIME_DURATION << "s\n\n";
	Log::print();

	return {
		vnode_index_to_share_cell_indexes,
		boudnary_oc_indexes, boundary_normals,
		periodic_boundary_oc_nc_index_pairs, periodic_boundary_normals,
		inner_face_oc_nc_index_pairs, inner_face_normals };
}


template<size_t space_dimension>
std::vector<size_t> Grid_Builder<space_dimension>::find_cell_indexes_have_these_vnodes(const std::unordered_map<size_t, std::set<size_t>>& vertex_node_index_to_cell_container_indexes, const std::vector<size_t>& face_node_indexes) {
	const auto start_node_index = face_node_indexes[0];
	const auto end_node_index = face_node_indexes[1];

	const auto& indexes_have_start_node = vertex_node_index_to_cell_container_indexes.at(start_node_index);
	const auto& indexes_have_end_node = vertex_node_index_to_cell_container_indexes.at(end_node_index);

	std::vector<size_t> cell_continaer_indexes_have_these_nodes;
	std::set_intersection(indexes_have_start_node.begin(), indexes_have_start_node.end(), indexes_have_end_node.begin(), indexes_have_end_node.end(), std::back_inserter(cell_continaer_indexes_have_these_nodes));

	return cell_continaer_indexes_have_these_nodes;
}