#pragma once
#include "Geometry.h"

#include <set>
#include <unordered_map>
#include <unordered_set>

struct Cell_Grid_Information
{
	std::vector<Physical_Domain_Vector> centers;
	std::vector<double> volumes;
	std::vector<std::array<double, PHYSICAL_DOMAIN_DIMENSION>> coordinate_projected_volumes;
};
struct Inner_Face_Grid_Information
{
	std::vector<double> areas;
	std::vector<Physical_Domain_Vector> normals;
	std::vector<std::pair<size_t, size_t>> owner_neighbor_container_indexes;
};
struct Boundary_Grid_Information
{
	std::vector<ElementType> boudnary_types;
	std::vector<double> areas;
	std::vector<Physical_Domain_Vector> normals;
	std::vector<size_t> owner_container_indexes;
};
struct Grid_Information
{
	Cell_Grid_Information cell_grid_information;
	Boundary_Grid_Information boundary_grid_information;
	Inner_Face_Grid_Information inner_face_grid_information;
};


class Grid_Data_Converter
{
public:
	static Grid_Information convert(const Grid_Data& grid_data);

private:
	Grid_Data_Converter(void) = delete;

	static std::vector<Physical_Domain_Vector> extract_by_index(const std::vector<Physical_Domain_Vector>& nodes, const std::vector<size_t>& indexes);
	static std::vector<size_t> find_cell_index_has_these_nodes(const std::unordered_map<size_t, std::set<size_t>>& vertex_node_index_to_cell_index, const std::vector<size_t>& face_node_indexes);
	static std::unordered_map<size_t, size_t> match_periodic_boudnary_index(const std::unordered_map<size_t, Geometry>& data_index_to_geometry, const ElementType element_type);
};


// template definition part
//template <typename R>
//Grid_Information Grid_Information_Builder<R>::build(const std::string& grid_file_path) {

