//#pragma once
//
//#include "Grid_Data_Processor.h"
//#include "Spatial_Discrete_Method.h"
//#include "Reconstruction_Method.h"
//
//
//template <size_t space_dimension>
//struct Cells_Info_FVM_Base
//{
//	std::vector<EuclideanVector<space_dimension>> centers;
//	std::vector<double> volumes;
//	std::vector<std::array<double, space_dimension>> coordinate_projected_volumes;
//};
//
//
//template <size_t space_dimension>
//struct Cells_Info_FVM_MLP
//{
//	Cells_Info_FVM_Base<space_dimension> cells_info_base;
//	std::vector<std::vector<size_t>> cells_vertex_node_indexes;
//	std::unordered_map<size_t, std::vector<size_t>> vertex_node_index_to_stencil;
//	std::vector<Dynamic_Matrix_> center_to_center_matrixes;
//	std::vector<Dynamic_Matrix_> center_to_vertex_matrixes;
//};
//
//
//template <size_t space_dimension>
//struct Inner_Faces_Info_FVM_Base
//{
//	std::vector<double> areas;
//	std::vector<EuclideanVector<space_dimension>> normals;
//	std::vector<std::pair<size_t, size_t>> owner_neighbor_container_indexes;
//};
//
//template <size_t space_dimension>
//struct Inner_Faces_Info_FVM_MLP
//{
//	using Space_Vector = EuclideanVector<space_dimension>;
//
//	Inner_Faces_Info_FVM_Base inner_faces_info_base;
//	std::vector<std::pair<Space_Vector, Space_Vector>> cell_to_face_vectors_o_n;
//};
//
//
//template <size_t space_dimension>
//struct Boundaries_Info_FVM
//{
//	std::vector<ElementType> types;
//	std::vector<double> areas;
//	std::vector<EuclideanVector<space_dimension>> normals;
//	std::vector<size_t> owner_container_indexes;
//};
//
//
//template <typename Spatial_Discrete_Method, typename Reconstruction_Method, size_t space_dimension>
//class Grid_Info_Extractor;
//
//
//template <size_t space_dimension>
//class Grid_Info_Extractor<FVM, Constant_Reconstruction, space_dimension>
//{
//public:
//	struct Grid_Infos
//	{
//		Cells_Info_FVM_Base<space_dimension> cells_information;
//		Boundaries_Info_FVM<space_dimension> boundaries_information;
//		Inner_Faces_Info_FVM_Base<space_dimension> inner_faces_information;
//	};
//
//	static Grid_Infos extract(Processed_Grid_Data<space_dimension>&& processed_grid_data);
//
//private:
//	Grid_Info_Extractor(void) = delete;
//};
//
//
//template <typename Gradient_Method, size_t space_dimension>
//class Grid_Info_Extractor<FVM, MLP_u1<Gradient_Method>, space_dimension>
//{
//public:
//	struct Grid_Infos
//	{
//		Cells_Info_FVM_MLP<space_dimension> cells_information;
//		Boundaries_Info_FVM<space_dimension> boundaries_information;
//		Inner_Faces_Info_FVM_MLP<space_dimension> inner_faces_information;
//	};
//
//	static Grid_Infos extract(Processed_Grid_Data<space_dimension>&& processed_grid_data);
//
//private:
//	Grid_Info_Extractor(void) = delete;
//};
//
//
////template definition part
//template <size_t space_dimension>
//Grid_Info_Extractor<FVM, Constant_Reconstruction, space_dimension>::Grid_Infos Grid_Info_Extractor<FVM, Constant_Reconstruction, space_dimension>::extract(Processed_Grid_Data<space_dimension>&& processed_grid_data) {
//	//// build cell infos
//	//const auto& cell_geometries = processed_grid_data.cell_geometries;
//	//const auto num_cell = cell_geometries.size();
//
//	//Cells_Info_FVM_Base<space_dimension> cells_info;
//	//cells_info.centers.reserve(num_cell);
//	//cells_info.coordinate_projected_volumes.reserve(num_cell);
//	//cells_info.volumes.reserve(num_cell);
//	//for (const auto& geometry : cell_geometries) {
//	//	cells_info.centers.push_back(geometry.center_node());
//	//	cells_info.coordinate_projected_volumes.push_back(geometry.coordinate_projected_volume());
//	//	cells_info.volumes.push_back(geometry.volume());
//	//}
//
//	// build boundary infos
//	const auto& boundary_geometries = processed_grid_data.boundary_geometries;
//	const auto num_boundary = boundary_geometries.size();
//
//	Boundaries_Info_FVM<space_dimension> boundaries_info;
//	boundaries_info.types = std::move(processed_grid_data.boundary_types);
//	boundaries_info.owner_container_indexes = std::move(processed_grid_data.boudnary_owner_container_indexes);
//	boundaries_info.areas.reserve(num_boundary);
//	boundaries_info.normals.reserve(num_boundary);
//	for (size_t i = 0; i < num_boundary; ++i) {
//		const auto& geometry = boundary_geometries[i];
//		const auto owner_cell_container_index = boundaries_info.owner_container_indexes[i];
//		const auto& owner_cell_center = cells_info.centers[owner_cell_container_index];
//
//		boundaries_info.areas.push_back(geometry.volume());
//		boundaries_info.normals.push_back(geometry.normal_vector(owner_cell_center));
//	}
//
//	// build inner face infos
//	// periodic boundary can be seen as inner face
//	const auto num_periodic_boundary_pair = processed_grid_data.periodic_boundary_owner_side_geometries.size();
//	const auto num_inner_face = processed_grid_data.inner_face_geometries.size();
//	const auto num_total_inner_face = num_periodic_boundary_pair + num_inner_face;
//
//	Inner_Faces_Info_FVM_Base<space_dimension> inner_faces_info;
//	inner_faces_info.areas.reserve(num_total_inner_face);
//	inner_faces_info.normals.reserve(num_total_inner_face);
//	inner_faces_info.owner_neighbor_container_indexes.reserve(num_total_inner_face);
//
//	for (size_t i = 0; i < num_periodic_boundary_pair; ++i) {
//		const auto& owner_side_geometry = processed_grid_data.periodic_boundary_owner_side_geometries[i];
//		const auto& [container_index_o, container_index_n] = processed_grid_data.periodic_boundary_owner_neighbor_container_indexes[i];
//		const auto owner_cell_center = cells_info.centers.at(container_index_o);
//
//		inner_faces_info.areas.push_back(owner_side_geometry.volume());
//		inner_faces_info.owner_neighbor_container_indexes.push_back({ container_index_o,container_index_n });
//		inner_faces_info.normals.push_back(owner_side_geometry.normal_vector(owner_cell_center));
//	}
//
//	for (size_t i = 0; i < num_inner_face; ++i) {
//		const auto& inner_face_geometry = processed_grid_data.inner_face_geometries[i];
//		const auto& [container_index_o, container_index_n] = processed_grid_data.inner_face_owner_neighbor_container_indexes[i];
//		const auto owner_cell_center = cells_info.centers.at(container_index_o);
//
//		inner_faces_info.areas.push_back(inner_face_geometry.volume());
//		inner_faces_info.owner_neighbor_container_indexes.push_back({ container_index_o,container_index_n });
//		inner_faces_info.normals.push_back(inner_face_geometry.normal_vector(owner_cell_center));
//	}
//
//	return { cells_info, boundaries_info, inner_faces_info };
//}
//
//
//template <typename Gradient_Method, size_t space_dimension>
//Grid_Info_Extractor<FVM, MLP_u1<Gradient_Method>, space_dimension>::Grid_Infos Grid_Info_Extractor<FVM, MLP_u1<Gradient_Method>, space_dimension>::extract(Processed_Grid_Data<space_dimension>&& processed_grid_data) {
//	// build cell infos
//	const auto& cell_geometries = processed_grid_data.cell_geometries;
//	const auto num_cell = cell_geometries.size();
//
//	Cells_Info_FVM_Base<space_dimension> cells_info;
//	cells_info.centers.reserve(num_cell);
//	cells_info.coordinate_projected_volumes.reserve(num_cell);
//	cells_info.volumes.reserve(num_cell);
//	for (const auto& geometry : cell_geometries) {
//		cells_info.centers.push_back(geometry.center_node());
//		cells_info.coordinate_projected_volumes.push_back(geometry.coordinate_projected_volume());
//		cells_info.volumes.push_back(geometry.volume());
//	}
//
//	// build boundary infos
//	const auto& boundary_geometries = processed_grid_data.boundary_geometries;
//	const auto num_boundary = boundary_geometries.size();
//
//	Boundaries_Info_FVM<space_dimension> boundaries_info;
//	boundaries_info.types = std::move(processed_grid_data.boundary_types);
//	boundaries_info.owner_container_indexes = std::move(processed_grid_data.boudnary_owner_container_indexes);
//	boundaries_info.areas.reserve(num_boundary);
//	boundaries_info.normals.reserve(num_boundary);
//	for (size_t i = 0; i < num_boundary; ++i) {
//		const auto& geometry = boundary_geometries[i];
//		const auto owner_cell_container_index = boundaries_info.owner_container_indexes[i];
//		const auto& owner_cell_center = cells_info.centers[owner_cell_container_index];
//
//		boundaries_info.areas.push_back(geometry.volume());
//		boundaries_info.normals.push_back(geometry.normal_vector(owner_cell_center));
//	}
//
//	// build inner face infos
//	// periodic boundary can be seen as inner face
//	const auto num_periodic_boundary_pair = processed_grid_data.periodic_boundary_owner_side_geometries.size();
//	const auto num_inner_face = processed_grid_data.inner_face_geometries.size();
//	const auto num_total_inner_face = num_periodic_boundary_pair + num_inner_face;
//
//	Inner_Faces_Info_FVM_Base<space_dimension> inner_faces_info;
//	inner_faces_info.areas.reserve(num_total_inner_face);
//	inner_faces_info.normals.reserve(num_total_inner_face);
//	inner_faces_info.owner_neighbor_container_indexes.reserve(num_total_inner_face);
//
//	for (size_t i = 0; i < num_periodic_boundary_pair; ++i) {
//		const auto& owner_side_geometry = processed_grid_data.periodic_boundary_owner_side_geometries[i];
//		const auto& [container_index_o, container_index_n] = processed_grid_data.periodic_boundary_owner_neighbor_container_indexes[i];
//		const auto owner_cell_center = cells_info.centers.at(container_index_o);
//
//		inner_faces_info.areas.push_back(owner_side_geometry.volume());
//		inner_faces_info.owner_neighbor_container_indexes.push_back({ container_index_o,container_index_n });
//		inner_faces_info.normals.push_back(owner_side_geometry.normal_vector(owner_cell_center));
//	}
//
//	for (size_t i = 0; i < num_inner_face; ++i) {
//		const auto& inner_face_geometry = processed_grid_data.inner_face_geometries[i];
//		const auto& [container_index_o, container_index_n] = processed_grid_data.inner_face_owner_neighbor_container_indexes[i];
//		const auto owner_cell_center = cells_info.centers.at(container_index_o);
//
//		inner_faces_info.areas.push_back(inner_face_geometry.volume());
//		inner_faces_info.owner_neighbor_container_indexes.push_back({ container_index_o,container_index_n });
//		inner_faces_info.normals.push_back(inner_face_geometry.normal_vector(owner_cell_center));
//	}
//
//	return { cells_info, boundaries_info, inner_faces_info };
//}
//
//
