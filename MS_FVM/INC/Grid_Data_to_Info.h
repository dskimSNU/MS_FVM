#pragma once

#include "Grid_File_to_Data.h"
#include "Post.h"

#include <set>
#include <unordered_map>
#include <unordered_set>


template<size_t dim>
class Grid_Data_to_Info
{
	using Space_Vector = EuclideanVector<dim>;

public:
	struct Cell_Info
	{
		std::vector<Space_Vector> centers;
		std::vector<double> volumes;
		std::vector<std::array<double, dim>> coordinate_projected_volumes;
	};


	struct Inner_Face_Info
	{
		std::vector<double> areas;
		std::vector<Space_Vector> normals;
		std::vector<std::pair<size_t, size_t>> owner_neighbor_container_indexes;
	};


	struct Boundary_Info
	{
		std::vector<ElementType> boudnary_types;
		std::vector<double> areas;
		std::vector<Space_Vector> normals;
		std::vector<size_t> owner_container_indexes;
	};


	struct Grid_Info
	{
		Cell_Info cell_grid_information;
		Boundary_Info boundary_grid_information;
		Inner_Face_Info inner_face_grid_information;
	};

	static Grid_Info convert(Grid_Data<dim>&& grid_data);

private:
	Grid_Data_to_Info(void) = delete;

	static std::vector<Space_Vector> extract_by_index(const std::vector<Space_Vector>& nodes, const std::vector<size_t>& indexes);
	static std::vector<size_t> find_cell_index_has_these_nodes(const std::unordered_map<size_t, std::set<size_t>>& vertex_node_index_to_cell_index, const std::vector<size_t>& face_node_indexes);
	static std::unordered_map<size_t, size_t> match_periodic_boudnary_index(const std::unordered_map<size_t, Geometry<dim>>& data_index_to_geometry, const ElementType element_type);
};


//template definition part
template <size_t dim>
Grid_Data_to_Info<dim>::Grid_Info Grid_Data_to_Info<dim>::convert(Grid_Data<dim>&& grid_data) {
	auto& [node_datas, cell_datas, boundary_datas, periodic_boundary_datas] = grid_data;

	// make 'vertex node index to cell index' to find owner/neighbor cell index
	const auto num_node = node_datas.size();
	std::unordered_map<size_t, std::set<size_t>> vertex_node_index_to_cell_index;
	vertex_node_index_to_cell_index.reserve(num_node);

	// make 'cell index to container index' to find owner/neighbor cell container index
	const auto num_cell = cell_datas.size();
	std::unordered_map<size_t, size_t> cell_index_to_container_index;
	cell_index_to_container_index.reserve(num_cell);

	// make 'cell_geometry'
	std::vector<Geometry<dim>> cell_geometries;
	cell_geometries.reserve(num_cell);

	// make 'faces_node_indexes' for construct inner face geometry datas
	std::set<Geometry<dim>> face_geometries;	//이걸 set에 넣으려고하니 Geometry 대소비교가 가능해야되고 그럴려니 EuclideanVector가 대소비교가 가능해야 된다.. 하 이게 뭐지


	// build cell_grid_info
	Cell_Info cell_grid_info;
	cell_grid_info.centers.resize(num_cell);
	cell_grid_info.volumes.resize(num_cell);
	cell_grid_info.coordinate_projected_volumes.resize(num_cell);
	for (size_t i = 0; i < num_cell; ++i) {
		auto& [index, figure, figure_order, type, node_indexes] = cell_datas[i];

		auto nodes = extract_by_index(node_datas, node_indexes);
		Geometry<dim> geometry(figure, figure_order, std::move(nodes), std::move(node_indexes));

		//
		cell_grid_info.centers[i] = geometry.center_node();
		cell_grid_info.volumes[i] = geometry.volume();
		cell_grid_info.coordinate_projected_volumes[i] = geometry.coordinate_projected_volume();

		//
		const auto vertex_node_indexes = geometry.vertex_node_indexes();
		for (const auto& vertex_node_index : vertex_node_indexes) {
			if (vertex_node_index_to_cell_index.find(vertex_node_index) == vertex_node_index_to_cell_index.end())
				vertex_node_index_to_cell_index.emplace(vertex_node_index, std::set<size_t>());

			vertex_node_index_to_cell_index.at(vertex_node_index).insert(index);
		}

		//
		auto cell_face_geometries = geometry.face_geometries();
		face_geometries.insert(std::make_move_iterator(cell_face_geometries.begin()), std::make_move_iterator(cell_face_geometries.end()));

		cell_index_to_container_index.emplace(index, i);
		cell_geometries.push_back(std::move(geometry));
	}

	
	// post grid
	Post::grid(cell_geometries);
	//std::vector<Space_Vector> global_post_nodes;
	//std::vector<std::vector<size_t>> global_connectivities;
	//for (const auto& geometry : cell_geometries) {
	//	auto connectivities = geometry.local_connectivities();
	//	auto post_nodes = geometry.vertex_nodes();	//post order를 geometry에게 어떻게 전달해줄 것인가?
	//	const auto num_post_node = post_nodes.size();

	//	Post::convert_to_global_connectivities(connectivities, num_post_node);

	//	global_post_nodes.insert(global_post_nodes.end(), std::make_move_iterator(post_nodes.begin()), std::make_move_iterator(post_nodes.end()));
	//	global_connectivities.insert(global_connectivities.end(), std::make_move_iterator(connectivities.begin()), std::make_move_iterator(connectivities.end()));
	//}
	//Post::write_grid_file(global_post_nodes, global_connectivities);


	// build boundary_grid_info
	Boundary_Info boundary_grid_info;
	if (!boundary_datas.empty()) {
		const auto num_boundary = boundary_datas.size();
		boundary_grid_info.boudnary_types.resize(num_boundary);
		boundary_grid_info.areas.resize(num_boundary);
		boundary_grid_info.normals.resize(num_boundary);
		boundary_grid_info.owner_container_indexes.resize(num_boundary);

		for (size_t i = 0; i < num_boundary; ++i) {
			auto& [index, figure, figure_order, type, node_indexes] = boundary_datas[i];

			auto nodes = extract_by_index(node_datas, node_indexes);
			Geometry geometry(figure, figure_order, std::move(nodes), std::move(node_indexes));

			//find owner cell container index
			const auto cell_indexes = find_cell_index_has_these_nodes(vertex_node_index_to_cell_index, node_indexes);
			dynamic_require(cell_indexes.size() == 1, "boundary should have unique owner cell");
			const auto owner_cell_index = cell_indexes.front();
			const auto owner_cell_container_index = cell_index_to_container_index.at(owner_cell_index);

			const auto& owner_cell_geometry = cell_geometries.at(owner_cell_container_index);

			//
			boundary_grid_info.boudnary_types[i] = type;
			boundary_grid_info.areas[i] = geometry.volume();
			boundary_grid_info.normals[i] = geometry.normal_vector(owner_cell_geometry.center_node());
			boundary_grid_info.owner_container_indexes[i] = owner_cell_container_index;

			//erase constructed face geometries
			const auto result = face_geometries.erase(geometry);
			dynamic_require(result == 1, "faces_geometries should have every boundary geometry");
		}
	}

	// build inner face grid info
	// periodic boundary can be seen as inner face
	Inner_Face_Info inner_face_grid_info;
	if (!periodic_boundary_datas.empty()) {
		const auto num_periodic_boundary = periodic_boundary_datas.size();
		const auto num_inner_face = static_cast<size_t>(num_periodic_boundary * 0.5);
		inner_face_grid_info.areas.reserve(num_periodic_boundary);
		inner_face_grid_info.normals.reserve(num_periodic_boundary);
		inner_face_grid_info.owner_neighbor_container_indexes.reserve(num_periodic_boundary);

		//build data index to geometry
		std::unordered_map<size_t, Geometry<dim>> data_index_to_x_periodic_geoemty, data_index_to_y_periodic_geoemty;
		for (size_t i = 0; i < num_periodic_boundary; ++i) {
			auto& [index, figure, figure_order, type, node_indexes] = periodic_boundary_datas[i];

			auto nodes = extract_by_index(node_datas, node_indexes);
			Geometry geometry(figure, figure_order, std::move(nodes), std::move(node_indexes));

			switch (type) {
			case ElementType::x_periodic:
				data_index_to_x_periodic_geoemty.emplace(i, std::move(geometry));
				break;
			case ElementType::y_periodic:
				data_index_to_y_periodic_geoemty.emplace(i, std::move(geometry));
				break;
			default:
				throw std::runtime_error("not supproted element type");
				break;
			}
		}

		std::unordered_map<size_t, size_t> data_index_to_matched_index;
		data_index_to_matched_index.reserve(num_inner_face);
		data_index_to_matched_index.merge(match_periodic_boudnary_index(data_index_to_x_periodic_geoemty, ElementType::x_periodic));
		data_index_to_matched_index.merge(match_periodic_boudnary_index(data_index_to_y_periodic_geoemty, ElementType::y_periodic));

		std::unordered_map<size_t, Geometry<dim>> data_index_to_geometry;
		data_index_to_geometry.reserve(num_periodic_boundary);
		data_index_to_geometry.merge(std::move(data_index_to_x_periodic_geoemty));
		data_index_to_geometry.merge(std::move(data_index_to_y_periodic_geoemty));

		//i : one of periodic face, j : matched periodic face of i
		for (const auto& [i_data_index, j_data_index] : data_index_to_matched_index) {
			const auto& owner_side_geometry = data_index_to_geometry.at(i_data_index);
			const auto& neighbor_side_geometry = data_index_to_geometry.at(j_data_index);

			//find owner/neighbor cell container indexes
			//we will designate i as owner cell side face		
			const auto cell_index_set_have_i = find_cell_index_has_these_nodes(vertex_node_index_to_cell_index, owner_side_geometry.vertex_node_indexes());
			const auto cell_index_set_have_j = find_cell_index_has_these_nodes(vertex_node_index_to_cell_index, neighbor_side_geometry.vertex_node_indexes());
			dynamic_require(cell_index_set_have_i.size() == 1, "periodic boundary should have unique owner cell");
			dynamic_require(cell_index_set_have_j.size() == 1, "periodic boundary should have unique neighbor cell");

			const auto owner_cell_index = cell_index_set_have_i.front();
			const auto neighbor_cell_index = cell_index_set_have_j.front();

			const auto owner_cell_container_index = cell_index_to_container_index.at(owner_cell_index);
			const auto neighbor_cell_container_index = cell_index_to_container_index.at(neighbor_cell_index);

			const auto& owner_cell_geomtry = cell_geometries.at(owner_cell_container_index);

			//			
			inner_face_grid_info.areas.push_back(owner_side_geometry.volume());
			inner_face_grid_info.normals.push_back(owner_side_geometry.normal_vector(owner_cell_geomtry.center_node()));
			inner_face_grid_info.owner_neighbor_container_indexes.push_back({ owner_cell_container_index,neighbor_cell_container_index });


			//erase constructed face geometries
			const auto i_result = face_geometries.erase(owner_side_geometry);
			const auto j_result = face_geometries.erase(neighbor_side_geometry);
			dynamic_require(i_result == 1, "face_geometries should have periodic boundary face");
			dynamic_require(j_result == 1, "face_geometries should have periodic boundary face");
		}
	}

	//build inner face
	for (const auto& face_geometry : face_geometries) {
		// need to know face geometry consisting node indexes
		const auto cell_indexes = find_cell_index_has_these_nodes(vertex_node_index_to_cell_index, face_geometry.vertex_node_indexes());
		dynamic_require(cell_indexes.size() == 2, "inner face should have owner neighbor cell");

		const auto owner_cell_index = cell_indexes[0];
		const auto neighbor_cell_index = cell_indexes[1];

		const auto owner_cell_container_indexes = cell_index_to_container_index.at(owner_cell_index);
		const auto neighbor_cell_container_indexes = cell_index_to_container_index.at(neighbor_cell_index);

		const auto& owner_cell_geometry = cell_geometries.at(owner_cell_container_indexes);

		inner_face_grid_info.areas.push_back(face_geometry.volume());
		inner_face_grid_info.normals.push_back(face_geometry.normal_vector(owner_cell_geometry.center_node()));
		inner_face_grid_info.owner_neighbor_container_indexes.push_back({ owner_cell_container_indexes,neighbor_cell_container_indexes });
	}

	return { cell_grid_info, boundary_grid_info, inner_face_grid_info };
}

template <size_t dim>
std::vector<EuclideanVector<dim>> Grid_Data_to_Info<dim>::extract_by_index(const std::vector<Space_Vector>& nodes, const std::vector<size_t>& indexes) {
	const auto num_extracted_node = indexes.size();
	std::vector<Space_Vector> extracted_nodes(num_extracted_node);
	for (size_t i = 0; i < num_extracted_node; ++i)
		extracted_nodes[i] = nodes[indexes[i] - 1]; // index start with 1

	return extracted_nodes;
}

template <size_t dim>
std::vector<size_t> Grid_Data_to_Info<dim>::find_cell_index_has_these_nodes(const std::unordered_map<size_t, std::set<size_t>>& vertex_node_index_to_cell_index, const std::vector<size_t>& face_node_indexes) {
	const auto start_node_index = face_node_indexes[0];
	const auto end_node_index = face_node_indexes[1];

	const auto& cell_index_set_have_start_node = vertex_node_index_to_cell_index.at(start_node_index);
	const auto& cell_index_set_have_end_node = vertex_node_index_to_cell_index.at(end_node_index);

	std::vector<size_t> cell_index_intersection;
	std::set_intersection(cell_index_set_have_start_node.begin(), cell_index_set_have_start_node.end(), cell_index_set_have_end_node.begin(), cell_index_set_have_end_node.end(), std::back_inserter(cell_index_intersection));

	return cell_index_intersection;
}

template <size_t dim>
std::unordered_map<size_t, size_t> Grid_Data_to_Info<dim>::match_periodic_boudnary_index(const std::unordered_map<size_t, Geometry<dim>>& data_index_to_geometry, const ElementType type) {
	const auto num_periodic_face = data_index_to_geometry.size();

	size_t axis_tag = 0;
	switch (type)
	{
	case ElementType::x_periodic:
		axis_tag = 0; break;
	case ElementType::y_periodic:
		axis_tag = 1; break;
	default:
		throw std::runtime_error("wrong element type");
		break;
	}

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

