#pragma once
#include "Grid_File_Convertor.h"
#include "Post.h"
#include "Profiler.h"

#include <set>
#include <unordered_map>
#include <unordered_set>


template <size_t space_dimension>
struct Processed_Grid_Data
{
	std::vector<Geometry<space_dimension>> cell_geometries;
	std::unordered_map<size_t, std::set<size_t>> vertex_node_index_to_cell_container_indexes;

	std::vector<ElementType> boundary_types;
	std::vector<Geometry<space_dimension>> boundary_geometries; // area, normal
	std::vector<size_t> boudnary_owner_container_indexes;

	std::vector<Geometry<space_dimension>> periodic_boundary_owner_side_geometries; //area, normal
	std::vector<Geometry<space_dimension>> periodic_boundary_neighbor_side_geometries;
	std::vector<std::pair<size_t, size_t>> periodic_boundary_owner_neighbor_container_indexes;

	std::vector<Geometry<space_dimension>> inner_face_geometries; //area, normal  
	std::vector<std::pair<size_t, size_t>> inner_face_owner_neighbor_container_indexes;
};


template<size_t space_dimension>
class Grid_Data_Processor
{
	using Space_Vector = EuclideanVector<space_dimension>;
public:
	static Processed_Grid_Data<space_dimension> process(Grid_Raw_Data<space_dimension>&& grid_data);

private:
	static void process_cell_data(Processed_Grid_Data<space_dimension>& process_grid_data, Grid_Raw_Data<space_dimension>&& grid_data);
	static void process_boudnary_data(Processed_Grid_Data<space_dimension>& process_grid_data, Grid_Raw_Data<space_dimension>&& grid_data);
	static void process_periodic_boundary_data(Processed_Grid_Data<space_dimension>& process_grid_data, Grid_Raw_Data<space_dimension>&& grid_data);
	static void process_inner_face_data(Processed_Grid_Data<space_dimension>& process_grid_data);


	static std::vector<Space_Vector> extract_by_index(const std::vector<Space_Vector>& nodes, const std::vector<size_t>& indexes);
	static std::vector<size_t> find_cell_container_indexes_have_these_nodes(const std::unordered_map<size_t, std::set<size_t>>& vertex_node_index_to_cell_container_indexes, const std::vector<size_t>& face_node_indexes);
	static std::unordered_map<size_t, size_t> match_periodic_boudnary_index(const std::unordered_map<size_t, Geometry<space_dimension>>& data_index_to_geometry, const ElementType element_type);
};


//template definition part
template <size_t space_dimension>
Processed_Grid_Data<space_dimension> Grid_Data_Processor<space_dimension>::process(Grid_Raw_Data<space_dimension>&& grid_data) {
	std::cout << std::left;
	std::cout << "============================================================\n";
	std::cout << "\t Process grid data\n";
	std::cout << "============================================================\n";
	SET_TIME_POINT;

	Processed_Grid_Data<space_dimension> processed_grid_data;
	process_cell_data(processed_grid_data, std::move(grid_data));
	process_boudnary_data(processed_grid_data, std::move(grid_data));
	process_periodic_boundary_data(processed_grid_data, std::move(grid_data));
	process_inner_face_data(processed_grid_data);

	std::cout << "============================================================\n";
	std::cout << "\t Total ellapsed time: " << std::setw(10) << GET_TIME_DURATION << "s\n";
	std::cout << "============================================================\n\n\n";

	// post grid
	Post::grid(processed_grid_data.cell_geometries);
	return processed_grid_data;
}


template <size_t space_dimension>
void Grid_Data_Processor<space_dimension>::process_cell_data(Processed_Grid_Data<space_dimension>& processed_grid_data, Grid_Raw_Data<space_dimension>&& grid_data) {
	SET_TIME_POINT;
	auto& [node_datas, cell_datas, boundary_datas, periodic_boundary_datas] = grid_data;

	auto& cell_geometries = processed_grid_data.cell_geometries;
	auto& vertex_node_index_to_cell_container_indexes = processed_grid_data.vertex_node_index_to_cell_container_indexes;

	const auto num_node = node_datas.size();
	const auto num_cell = cell_datas.size();
	cell_geometries.reserve(num_cell);
	vertex_node_index_to_cell_container_indexes.reserve(num_node);

	for (size_t i = 0; i < num_cell; ++i) {
		auto& [index, figure, figure_order, type, node_indexes] = cell_datas[i];

		auto nodes = extract_by_index(node_datas, node_indexes);
		Geometry<space_dimension> geometry(figure, figure_order, std::move(nodes), std::move(node_indexes));

		const auto vertex_node_indexes = geometry.vertex_node_indexes();
		for (const auto& vertex_node_index : vertex_node_indexes) {
			if (vertex_node_index_to_cell_container_indexes.find(vertex_node_index) == vertex_node_index_to_cell_container_indexes.end())
				vertex_node_index_to_cell_container_indexes.emplace(vertex_node_index, std::set<size_t>());

			vertex_node_index_to_cell_container_indexes.at(vertex_node_index).insert(i);
		}

		cell_geometries.push_back(std::move(geometry));
	}
	std::cout << "process " << std::setw(5) << num_cell << " cells \t\t\t" << "ellapsed " << std::setw(10) << GET_TIME_DURATION << "s\n";
}


template <size_t space_dimension>
void Grid_Data_Processor<space_dimension>::process_boudnary_data(Processed_Grid_Data<space_dimension>& processed_grid_data, Grid_Raw_Data<space_dimension>&& grid_data) {
	SET_TIME_POINT;
	auto& [node_datas, cell_datas, boundary_datas, periodic_boundary_datas] = grid_data;

	if (boundary_datas.empty())
		return;

	//processed cell data
	auto& vertex_node_index_to_cell_container_indexes = processed_grid_data.vertex_node_index_to_cell_container_indexes;

	//to be processed boundary data
	const auto num_boundary = boundary_datas.size();

	auto& types = processed_grid_data.boundary_types;
	auto& owner_container_indexes = processed_grid_data.boudnary_owner_container_indexes;
	auto& geometries = processed_grid_data.boundary_geometries;

	types.resize(num_boundary);
	owner_container_indexes.resize(num_boundary);
	geometries.reserve(num_boundary);

	for (size_t i = 0; i < num_boundary; ++i) {
		auto& [index, figure, figure_order, type, node_indexes] = boundary_datas[i];

		auto nodes = extract_by_index(node_datas, node_indexes);
		Geometry geometry(figure, figure_order, std::move(nodes), std::move(node_indexes));

		//find owner cell container index
		const auto cell_container_indexes = find_cell_container_indexes_have_these_nodes(vertex_node_index_to_cell_container_indexes, node_indexes);
		dynamic_require(cell_container_indexes.size() == 1, "boundary should have unique owner cell");
		const auto owner_cell_container_index = cell_container_indexes.front();

		types[i] = type;
		owner_container_indexes[i] = owner_cell_container_index;
		geometries.push_back(std::move(geometry));
	}
	std::cout << "process " << std::setw(5) << num_boundary << " boundaries \t\t" << "ellapsed " << std::setw(10) << GET_TIME_DURATION << "s\n";
}


template <size_t space_dimension>
void Grid_Data_Processor<space_dimension>::process_periodic_boundary_data(Processed_Grid_Data<space_dimension>& processed_grid_data, Grid_Raw_Data<space_dimension>&& grid_data) {
	SET_TIME_POINT;
	auto& [node_datas, cell_datas, boundary_datas, periodic_boundary_datas] = grid_data;

	if (periodic_boundary_datas.empty())
		return;

	//processed cell data
	auto& vertex_node_index_to_cell_container_indexes = processed_grid_data.vertex_node_index_to_cell_container_indexes;

	//to be processed periodic boundary data
	const auto num_periodic_boundary = periodic_boundary_datas.size();
	const auto num_periodic_pair = static_cast<size_t>(num_periodic_boundary * 0.5);

	auto& owner_side_geometries = processed_grid_data.periodic_boundary_owner_side_geometries;
	auto& neighbor_side_geometries = processed_grid_data.periodic_boundary_neighbor_side_geometries;
	auto& cell_container_index_o_n = processed_grid_data.periodic_boundary_owner_neighbor_container_indexes;

	owner_side_geometries.reserve(num_periodic_boundary);
	neighbor_side_geometries.reserve(num_periodic_boundary);
	cell_container_index_o_n.reserve(num_periodic_boundary);

	//build data index to geometry
	std::unordered_map<size_t, Geometry<space_dimension>> data_index_to_x_periodic_geoemty, data_index_to_y_periodic_geoemty;
	for (size_t i = 0; i < num_periodic_boundary; ++i) {
		auto& [index, figure, figure_order, type, node_indexes] = periodic_boundary_datas[i];

		auto nodes = extract_by_index(node_datas, node_indexes);
		Geometry geometry(figure, figure_order, std::move(nodes), std::move(node_indexes));

		if (type == ElementType::x_periodic)
			data_index_to_x_periodic_geoemty.emplace(i, std::move(geometry));
		else if (type == ElementType::y_periodic)
			data_index_to_y_periodic_geoemty.emplace(i, std::move(geometry));
		else
			throw std::runtime_error("not supproted element type");
	}

	std::unordered_map<size_t, size_t> data_index_to_matched_index;
	data_index_to_matched_index.reserve(num_periodic_pair);
	data_index_to_matched_index.merge(match_periodic_boudnary_index(data_index_to_x_periodic_geoemty, ElementType::x_periodic));
	data_index_to_matched_index.merge(match_periodic_boudnary_index(data_index_to_y_periodic_geoemty, ElementType::y_periodic));

	std::unordered_map<size_t, Geometry<space_dimension>> data_index_to_geometry;
	data_index_to_geometry.reserve(num_periodic_boundary);
	data_index_to_geometry.merge(std::move(data_index_to_x_periodic_geoemty));
	data_index_to_geometry.merge(std::move(data_index_to_y_periodic_geoemty));

	//i : one of periodic face, j : matched periodic face of i
	for (const auto& [i_data_index, j_data_index] : data_index_to_matched_index) {
		auto& owner_side_geometry = data_index_to_geometry.at(i_data_index);
		auto& neighbor_side_geometry = data_index_to_geometry.at(j_data_index);

		//find owner/neighbor cell container indexes
		//we will designate i as owner cell side face		
		const auto container_indexes_set_have_i = find_cell_container_indexes_have_these_nodes(vertex_node_index_to_cell_container_indexes, owner_side_geometry.vertex_node_indexes());
		const auto container_indexes_set_have_j = find_cell_container_indexes_have_these_nodes(vertex_node_index_to_cell_container_indexes, neighbor_side_geometry.vertex_node_indexes());
		dynamic_require(container_indexes_set_have_i.size() == 1, "periodic boundary should have unique owner cell");
		dynamic_require(container_indexes_set_have_j.size() == 1, "periodic boundary should have unique neighbor cell");

		const auto cell_container_index_o = container_indexes_set_have_i.front();
		const auto cell_container_index_n = container_indexes_set_have_j.front();

		owner_side_geometries.push_back(std::move(owner_side_geometry));
		neighbor_side_geometries.push_back(std::move(neighbor_side_geometry));
		cell_container_index_o_n.push_back({ cell_container_index_o,cell_container_index_n });
	}

	// update vertex_node_index_to_cell_container_indexes
	for (size_t i = 0; i < num_periodic_pair; ++i) {
		const auto [cell_container_index_o, cell_container_index_n] = cell_container_index_o_n[i];
		const auto& owner_side_geometry = owner_side_geometries[i];
		const auto& neighbor_side_geometry = neighbor_side_geometries[i];

		const auto owner_side_vertex_node_indexes = owner_side_geometry.vertex_node_indexes();
		const auto neighbor_side_vertex_node_indexes = neighbor_side_geometry.vertex_node_indexes();

		for (const auto owner_side_vertex_node_index : owner_side_vertex_node_indexes) {
			auto& cell_container_indexes = vertex_node_index_to_cell_container_indexes.at(owner_side_vertex_node_index);
			cell_container_indexes.insert(cell_container_index_n);
		}
		for (const auto neighbor_side_vertex_node_index : neighbor_side_vertex_node_indexes) {
			auto& cell_container_indexes = vertex_node_index_to_cell_container_indexes.at(neighbor_side_vertex_node_index);
			cell_container_indexes.insert(cell_container_index_o);
		}
	}

	for (size_t i = 0; i < num_periodic_pair; ++i) {
		const auto [cell_container_index_o, cell_container_index_n] = cell_container_index_o_n[i];
		const auto& owner_side_geometry = owner_side_geometries[i];
		const auto& neighbor_side_geometry = neighbor_side_geometries[i];


	}


	std::cout << "process " << std::setw(5) << num_periodic_boundary << " periodic boundaries \t" << "ellapsed " << std::setw(10) << GET_TIME_DURATION << "s\n";
}


template <size_t space_dimension>
void Grid_Data_Processor<space_dimension>::process_inner_face_data(Processed_Grid_Data<space_dimension>& processed_grid_data) {
	SET_TIME_POINT;
	//processed cell data
	const auto& cell_geometries = processed_grid_data.cell_geometries;
	const auto& vertex_node_index_to_cell_container_indexes = processed_grid_data.vertex_node_index_to_cell_container_indexes;

	//processed bdry data
	const auto& boundary_geometries = processed_grid_data.boundary_geometries;

	//processed pbdry data
	const auto& periodic_boundary_owner_side_geometries = processed_grid_data.periodic_boundary_owner_side_geometries;
	const auto& periodic_boundary_neighbor_side_geometries = processed_grid_data.periodic_boundary_neighbor_side_geometries;

	//construct all face geometry
	std::set<Geometry<space_dimension>> face_geometries;
	for (const auto& geometry : cell_geometries) {
		auto cell_face_geometries = geometry.face_geometries();
		face_geometries.insert(std::make_move_iterator(cell_face_geometries.begin()), std::make_move_iterator(cell_face_geometries.end()));
	}

	//erase constructed face geometry
	for (const auto& boundray_geometry : boundary_geometries) {
		const auto result = face_geometries.erase(boundray_geometry);
		dynamic_require(result == 1, "faces_geometries should have every boundary geometry");
	}
	for (const auto& periodic_boundary_geometry : periodic_boundary_owner_side_geometries) {
		const auto result = face_geometries.erase(periodic_boundary_geometry);
		dynamic_require(result == 1, "faces_geometries should have every periodic boundary geometry");
	}
	for (const auto& periodic_boundary_geometry : periodic_boundary_neighbor_side_geometries) {
		const auto result = face_geometries.erase(periodic_boundary_geometry);
		dynamic_require(result == 1, "faces_geometries should have every periodic boundary geometry");
	}

	//to be processed inner face data
	const auto num_inner_face = face_geometries.size();

	auto& inner_face_geometries = processed_grid_data.inner_face_geometries;
	auto& owner_neighbor_container_indexes = processed_grid_data.inner_face_owner_neighbor_container_indexes;

	inner_face_geometries.reserve(num_inner_face);
	owner_neighbor_container_indexes.reserve(num_inner_face);

	for (auto& face_geometry : face_geometries) {
		const auto cell_container_indexes = find_cell_container_indexes_have_these_nodes(vertex_node_index_to_cell_container_indexes, face_geometry.vertex_node_indexes());
		dynamic_require(cell_container_indexes.size() == 2, "inner face should have owner neighbor cell");

		const auto container_index_o = cell_container_indexes[0];
		const auto container_index_n = cell_container_indexes[1];

		inner_face_geometries.push_back(std::move(face_geometry));
		owner_neighbor_container_indexes.push_back({ container_index_o, container_index_n });
	}
	std::cout << "process " << std::setw(5) << face_geometries.size() << " inner faces \t\t" << "ellapsed " << std::setw(10) << GET_TIME_DURATION << "s\n";
}

template<size_t space_dimension>
std::vector<size_t> Grid_Data_Processor<space_dimension>::find_cell_container_indexes_have_these_nodes(const std::unordered_map<size_t, std::set<size_t>>& vertex_node_index_to_cell_container_indexes, const std::vector<size_t>& face_node_indexes) {
	const auto start_node_index = face_node_indexes[0];
	const auto end_node_index = face_node_indexes[1];

	const auto& indexes_have_start_node = vertex_node_index_to_cell_container_indexes.at(start_node_index);
	const auto& indexes_have_end_node = vertex_node_index_to_cell_container_indexes.at(end_node_index);

	std::vector<size_t> cell_continaer_indexes_have_these_nodes;
	std::set_intersection(indexes_have_start_node.begin(), indexes_have_start_node.end(), indexes_have_end_node.begin(), indexes_have_end_node.end(), std::back_inserter(cell_continaer_indexes_have_these_nodes));

	return cell_continaer_indexes_have_these_nodes;
}

template <size_t space_dimension>
std::vector<EuclideanVector<space_dimension>> Grid_Data_Processor<space_dimension>::extract_by_index(const std::vector<Space_Vector>& nodes, const std::vector<size_t>& indexes) {
	const auto num_extracted_node = indexes.size();
	std::vector<Space_Vector> extracted_nodes(num_extracted_node);
	for (size_t i = 0; i < num_extracted_node; ++i)
		extracted_nodes[i] = nodes[indexes[i] - 1]; // index start with 1

	return extracted_nodes;
}



template <size_t space_dimension>
std::unordered_map<size_t, size_t> Grid_Data_Processor<space_dimension>::match_periodic_boudnary_index(const std::unordered_map<size_t, Geometry<space_dimension>>& data_index_to_geometry, const ElementType type) {
	const auto num_periodic_face = data_index_to_geometry.size();

	size_t axis_tag = 0;
	if (type == ElementType::x_periodic)
		axis_tag = 0;
	else if (type == ElementType::y_periodic)
		axis_tag = 1;
	else
		throw std::runtime_error("wrong element type");

	if (num_periodic_face == 0)
		return std::unordered_map<size_t, size_t>();

	std::unordered_set<size_t> matched_index;
	matched_index.reserve(num_periodic_face);

	const auto num_set = static_cast<size_t>(0.5 * num_periodic_face);

	std::unordered_map<size_t, size_t> index_to_matched_index;
	index_to_matched_index.reserve(num_set);

	const auto start_iter = data_index_to_geometry.begin();
	const auto end_iter = data_index_to_geometry.end();
	for (auto i_iter = start_iter; i_iter != end_iter; ++i_iter) {
		const auto& [i_index, i_geometry] = *i_iter;
		if (matched_index.find(i_index) != matched_index.end())
			continue;

		for (auto j_iter = std::next(i_iter, 1); j_iter != end_iter; ++j_iter) {
			const auto& [j_index, j_geometry] = *j_iter;

			if (i_geometry.is_periodic_pair(j_geometry, axis_tag)) {
				index_to_matched_index.emplace(i_index, j_index);
				matched_index.emplace(i_index);
				matched_index.emplace(j_index);
				break;
			}
		}
	}

	dynamic_require(index_to_matched_index.size() == num_set, "some periodic boundaries matched yet");
	return index_to_matched_index;
}
