#include "../INC/Grid_Information_Builder.h"

Grid_Information Grid_Data_Converter::convert(const Grid_Data& grid_data) {
	auto [node_grid_datas, cell_grid_datas, boundary_grid_datas, periodic_boundary_grid_datas] = grid_data;

	// make 'vertex node index to cell index' to find owner/neighbor cell index
	const auto num_node = node_grid_datas.size();
	std::unordered_map<size_t, std::set<size_t>> vertex_node_index_to_cell_index;
	vertex_node_index_to_cell_index.reserve(num_node);

	// make 'cell index to container index' to find owner/neighbor cell container index
	const auto num_cell = cell_grid_datas.size();
	std::unordered_map<size_t, size_t> cell_index_to_container_index;
	cell_index_to_container_index.reserve(num_cell);

	// make 'cell_geometry'
	std::vector<Geometry> cell_geometries;
	cell_geometries.reserve(num_cell);

	// make 'faces_node_indexes' for construct inner face geometry datas
	std::set<Geometry> face_geometries;	//이걸 set에 넣으려고하니 Geometry 대소비교가 가능해야되고 그럴려니 EuclideanVector가 대소비교가 가능해야 된다.. 하 이게 뭐지


	// build cell_grid_info
	Cell_Grid_Information cell_grid_info;	
	cell_grid_info.centers.resize(num_cell);
	cell_grid_info.volumes.resize(num_cell);
	cell_grid_info.coordinate_projected_volumes.resize(num_cell);
	for (size_t i = 0; i < num_cell; ++i) {
		auto& [index, figure, figure_order, type, node_indexes] = cell_grid_datas[i];

		ReferenceGeometry ref_geometry(figure, figure_order);
		auto consisting_nodes = Grid_Data_Converter::extract_by_index(node_grid_datas, node_indexes);
		Geometry geometry(ref_geometry, std::move(consisting_nodes), std::move(node_indexes));
		
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


	// build boundary_grid_info
	Boundary_Grid_Information boundary_grid_info;
	if (!boundary_grid_datas.empty()) {
		const auto num_boundary = boundary_grid_datas.size();		
		boundary_grid_info.boudnary_types.resize(num_boundary);
		boundary_grid_info.areas.resize(num_boundary);
		boundary_grid_info.normals.resize(num_boundary);
		boundary_grid_info.owner_container_indexes.resize(num_boundary);

		for (size_t i = 0; i < num_boundary; ++i) {
			auto& [index, figure, figure_order, type, node_indexes] = boundary_grid_datas[i];
			
			ReferenceGeometry ref_geometry(figure, figure_order);
			auto consisting_nodes = Grid_Data_Converter::extract_by_index(node_grid_datas, node_indexes);
			Geometry geometry(ref_geometry, std::move(consisting_nodes), std::move(node_indexes));

			//find owner cell container index
			const auto cell_indexes = Grid_Data_Converter::find_cell_index_has_these_nodes(vertex_node_index_to_cell_index, node_indexes);
			dynamic_require(cell_indexes.size() == 1, std::to_string(index) + " boundary does not have unique owner cell");
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
			dynamic_require(result == 1, "faces_geometries is not complete");
		}
	}

	// build inner face grid info
	// periodic boundary can be seen as inner face
	Inner_Face_Grid_Information inner_face_grid_info;
	if (!periodic_boundary_grid_datas.empty()) {
		const auto num_periodic_boundary = periodic_boundary_grid_datas.size();
		const auto num_inner_face = static_cast<size_t>(num_periodic_boundary * 0.5);		
		inner_face_grid_info.areas.reserve(num_periodic_boundary);
		inner_face_grid_info.normals.reserve(num_periodic_boundary);
		inner_face_grid_info.owner_neighbor_container_indexes.reserve(num_periodic_boundary);

		//build data index to geometry
		std::unordered_map<size_t, Geometry> data_index_to_x_periodic_geoemty, data_index_to_y_periodic_geoemty;
		for (size_t i = 0; i < num_periodic_boundary; ++i) {
			auto& [index, figure, figure_order, type, node_indexes] = periodic_boundary_grid_datas[i];

			ReferenceGeometry ref_geometry(figure, figure_order);
			auto consisting_nodes = Grid_Data_Converter::extract_by_index(node_grid_datas, node_indexes);
			Geometry geometry(ref_geometry, std::move(consisting_nodes), std::move(node_indexes));

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
		data_index_to_matched_index.merge(Grid_Data_Converter::match_periodic_boudnary_index(data_index_to_x_periodic_geoemty, ElementType::x_periodic));
		data_index_to_matched_index.merge(Grid_Data_Converter::match_periodic_boudnary_index(data_index_to_y_periodic_geoemty, ElementType::y_periodic));

		std::unordered_map<size_t, Geometry> data_index_to_geometry;
		data_index_to_geometry.reserve(num_periodic_boundary);
		data_index_to_geometry.merge(std::move(data_index_to_x_periodic_geoemty));
		data_index_to_geometry.merge(std::move(data_index_to_y_periodic_geoemty));

		//i : one of periodic face, j : matched periodic face of i
		for (const auto& [i_data_index, j_data_index] : data_index_to_matched_index) {
			const auto& owner_side_geometry = data_index_to_geometry.at(i_data_index);
			const auto& neighbor_side_geometry = data_index_to_geometry.at(j_data_index);

			//find owner/neighbor cell container indexes
			//we will designate i as owner cell side face		
			const auto cell_index_set_have_i = Grid_Data_Converter::find_cell_index_has_these_nodes(vertex_node_index_to_cell_index, owner_side_geometry.vertex_node_indexes());
			const auto cell_index_set_have_j = Grid_Data_Converter::find_cell_index_has_these_nodes(vertex_node_index_to_cell_index, neighbor_side_geometry.vertex_node_indexes());
			dynamic_require(cell_index_set_have_i.size() == 1, " periodic boundary have not unique owner cell");
			dynamic_require(cell_index_set_have_j.size() == 1, " periodic boundary have not unique neighbor cell");

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
			dynamic_require(i_result == 1, "face_geometries is not complete");
			dynamic_require(j_result == 1, "face_geometries is not complete");
		}
	}

	//build inner face
	for (const auto& face_geometry : face_geometries) { 
		// need to know face geometry consisting node indexes
		const auto cell_indexes = Grid_Data_Converter::find_cell_index_has_these_nodes(vertex_node_index_to_cell_index, face_geometry.vertex_node_indexes());
		dynamic_require(cell_indexes.size() == 2, "inner face does not have unique owner neighbor cell");

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

std::vector<Physical_Domain_Vector> Grid_Data_Converter::extract_by_index(const std::vector<Physical_Domain_Vector>& nodes, const std::vector<size_t>& indexes) {
	const auto num_consisting_node = indexes.size();
	std::vector<Physical_Domain_Vector> consisting_nodes(num_consisting_node);
	for (size_t i = 0; i < num_consisting_node; ++i)
		consisting_nodes[i] = nodes[indexes[i] - 1]; // index start with 1

	return consisting_nodes;
}

std::vector<size_t> Grid_Data_Converter::find_cell_index_has_these_nodes(const std::unordered_map<size_t, std::set<size_t>>& vertex_node_index_to_cell_index, const std::vector<size_t>& face_node_indexes) {
	const auto start_node_index = face_node_indexes[0];
	const auto end_node_index = face_node_indexes[1];

	const auto& cell_index_set_have_start_node = vertex_node_index_to_cell_index.at(start_node_index);
	const auto& cell_index_set_have_end_node = vertex_node_index_to_cell_index.at(end_node_index);

	std::vector<size_t> cell_index_intersection;
	std::set_intersection(cell_index_set_have_start_node.begin(), cell_index_set_have_start_node.end(), cell_index_set_have_end_node.begin(), cell_index_set_have_end_node.end(), std::back_inserter(cell_index_intersection));

	return cell_index_intersection;
}

std::unordered_map<size_t, size_t> Grid_Data_Converter::match_periodic_boudnary_index(const std::unordered_map<size_t, Geometry>& data_index_to_geometry, const ElementType type){
	const auto num_periodic_face = data_index_to_geometry.size();

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

			if (i_geometry.is_periodic_pair(j_geometry,type)) {
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