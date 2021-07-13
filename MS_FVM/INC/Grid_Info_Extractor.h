#pragma once

#include "Grid_Data_Processor.h"
#include "Spatial_Discrete_Method.h"


template <typename SDM, size_t dim>
class Grid_Info_Extractor;


template <size_t dim>
class Grid_Info_Extractor<FVM, dim> 
{
	using Space_Vector = EuclideanVector<dim>;
public:
	struct Cell_Infos
	{
		std::vector<Space_Vector> centers;
		std::vector<double> volumes;
		std::vector<std::array<double, dim>> coordinate_projected_volumes;
	};
	struct Inner_Face_Infos
	{
		std::vector<double> areas;
		std::vector<Space_Vector> normals;
		std::vector<std::pair<size_t, size_t>> owner_neighbor_container_indexes;
	};
	struct Boundary_Infos
	{
		std::vector<ElementType> types;
		std::vector<double> areas;
		std::vector<Space_Vector> normals;
		std::vector<size_t> owner_container_indexes;
	};
	struct Grid_Infos
	{
		Cell_Infos cell_informations;
		Boundary_Infos boundary_informations;
		Inner_Face_Infos inner_face_informations;
	};

	static Grid_Infos extract(Processed_Grid_Data<dim>&& processed_grid_data);

private:
	Grid_Info_Extractor(void) = delete;
};


//template definition part
template <size_t dim>
Grid_Info_Extractor<FVM, dim>::Grid_Infos Grid_Info_Extractor<FVM, dim>::extract(Processed_Grid_Data<dim>&& processed_grid_data) {
	// build cell infos
	const auto& cell_geometries = processed_grid_data.cell_geometries;
	const auto num_cell = cell_geometries.size();

	Cell_Infos cell_infos;
	cell_infos.centers.reserve(num_cell);
	cell_infos.coordinate_projected_volumes.reserve(num_cell);
	cell_infos.volumes.reserve(num_cell);
	for (const auto& geometry : cell_geometries) {
		cell_infos.centers.push_back(geometry.center_node());
		cell_infos.coordinate_projected_volumes.push_back(geometry.coordinate_projected_volume());
		cell_infos.volumes.push_back(geometry.volume());
	}

	// build boundary infos
	const auto& boundary_geometries = processed_grid_data.boundary_geometries;
	const auto num_boundary = boundary_geometries.size();

	Boundary_Infos boundary_infos;
	boundary_infos.types = std::move(processed_grid_data.boundary_types);
	boundary_infos.owner_container_indexes = std::move(processed_grid_data.boudnary_owner_container_indexes);
	boundary_infos.areas.reserve(num_boundary);
	boundary_infos.normals.reserve(num_boundary);
	for (size_t i = 0; i < num_boundary; ++i) {
		const auto& geometry = boundary_geometries[i];
		const auto owner_cell_container_index = boundary_infos.owner_container_indexes[i];
		const auto& owner_cell_center = cell_infos.centers[owner_cell_container_index];

		boundary_infos.areas.push_back(geometry.volume());
		boundary_infos.normals.push_back(geometry.normal_vector(owner_cell_center));
	}

	// build inner face infos
	// periodic boundary can be seen as inner face
	const auto num_periodic_boundary_pair = processed_grid_data.periodic_boundary_owner_side_geometries.size();
	const auto num_inner_face = processed_grid_data.inner_face_geometries.size();
	const auto num_total_inner_face = num_periodic_boundary_pair + num_inner_face;

	Inner_Face_Infos inner_face_infos;
	inner_face_infos.areas.reserve(num_total_inner_face);
	inner_face_infos.normals.reserve(num_total_inner_face);
	inner_face_infos.owner_neighbor_container_indexes.reserve(num_total_inner_face);

	for (size_t i = 0; i < num_periodic_boundary_pair; ++i) {
		const auto& owner_side_geometry = processed_grid_data.periodic_boundary_owner_side_geometries[i];
		const auto& [container_index_o, container_index_n] = processed_grid_data.periodic_boundary_owner_neighbor_container_indexes[i];
		const auto owner_cell_center = cell_infos.centers.at(container_index_o);

		inner_face_infos.areas.push_back(owner_side_geometry.volume());
		inner_face_infos.owner_neighbor_container_indexes.push_back({ container_index_o,container_index_n });
		inner_face_infos.normals.push_back(owner_side_geometry.normal_vector(owner_cell_center));
	}

	for (size_t i = 0; i < num_inner_face; ++i) {
		const auto& inner_face_geometry = processed_grid_data.inner_face_geometries[i];
		const auto& [container_index_o, container_index_n] = processed_grid_data.inner_face_owner_neighbor_container_indexes[i];
		const auto owner_cell_center = cell_infos.centers.at(container_index_o);

		inner_face_infos.areas.push_back(inner_face_geometry.volume());
		inner_face_infos.owner_neighbor_container_indexes.push_back({ container_index_o,container_index_n });
		inner_face_infos.normals.push_back(inner_face_geometry.normal_vector(owner_cell_center));
	}

	return { cell_infos, boundary_infos, inner_face_infos };
}