#include "../INC/Grid_Information_Builder.h"

Grid_Information Grid_Information_Builder::build(const Grid_Data& grid_data) {
	auto [node_grid_datas, cell_grid_datas, boundary_grid_datas, periodic_boudnary_grid_datas] = grid_data;

	// build cell_grid_information
	const auto num_node = node_grid_datas.size();
	const auto num_cell = cell_grid_datas.size();

	std::vector<double> volumes(num_cell);
	std::vector<std::array<double, s_physical_domain_dimension>> coordinate_projected_volumes(num_cell);

	// make 'cell index to cell set index' to find owner/neighbor cells index
	std::unordered_map<size_t, size_t> cell_index_to_cells_index;
	cell_index_to_cells_index.reserve(num_cell);

	// make 'vertex node index to cell index' to find owner/neighbor cell index
	std::unordered_map<size_t, std::set<size_t>> vertex_node_index_to_cell_index;
	vertex_node_index_to_cell_index.reserve(num_node);

	for (size_t i = 0; i < num_cell; ++i) {
		const auto [index, figure, figure_order, type, node_indexes] = cell_grid_datas[i];

		ReferenceGeometry ref_geometry(figure, figure_order);
		auto consisting_nodes = Grid_Information_Builder::extract_by_index(node_grid_datas, node_indexes);
		Geometry geometry(ref_geometry, std::move(consisting_nodes));

		volumes[i] = geometry.calculate_volume();
		coordinate_projected_volumes[i] = geometry.calculate_coordinate_projected_volume();

		cell_index_to_cells_index.emplace(index, i);

		const auto vertex_node_indexes = geometry.calculate_vertex_node_indexes(node_indexes);
		for (const auto& vertex_node_index : vertex_node_indexes) {
			if (vertex_node_index_to_cell_index.find(vertex_node_index) == vertex_node_index_to_cell_index.end())
				vertex_node_index_to_cell_index.emplace(vertex_node_index, std::set<size_t>());

			vertex_node_index_to_cell_index.at(vertex_node_index).insert(index);
		}
	}

	return Grid_Information();
}

std::vector<Physical_Domain_Vector> Grid_Information_Builder::extract_by_index(const std::vector<Physical_Domain_Vector>& nodes, const std::vector<size_t>& indexes) {
	const auto num_consisting_node = indexes.size();
	std::vector<Physical_Domain_Vector> consisting_nodes(num_consisting_node);
	for (size_t i = 0; i < num_consisting_node; ++i)
		consisting_nodes[i] = nodes[indexes[i]];

	return consisting_nodes;
}