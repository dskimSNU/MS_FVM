#pragma once
#include "Geometry.h"

#include <set>
#include <unordered_map>


struct Cell_Grid_Information
{
	std::vector<double> volumes;
	std::vector<std::array<double, s_physical_domain_dimension>> coordinate_projected_volumes;
};
struct Inner_Face_Grid_Information
{
	std::vector<double> areas;
	std::vector<Physical_Domain_Vector> normals;
	std::vector<std::pair<size_t, size_t>> owner_neighbor_indexes;
};
struct Boundary_Face_Grid_Information
{
	std::vector<ElementType> boudnary_types;
	std::vector<double> areas;
	std::vector<Physical_Domain_Vector> normals;
	std::vector<size_t> owner_indexes;
};
struct Grid_Information
{
	Cell_Grid_Information cell_grid_information;
	Inner_Face_Grid_Information inner_face_grid_information;
	Boundary_Face_Grid_Information boundary_face_grid_information;
};


class Grid_Information_Builder
{
public:
	static Grid_Information build(const Grid_Data& grid_data);

private:
	static std::vector<Physical_Domain_Vector> extract_by_index(const std::vector<Physical_Domain_Vector>& nodes, const std::vector<size_t>& indexes);
};


// template definition part
//template <typename R>
//Grid_Information Grid_Information_Builder<R>::build(const std::string& grid_file_path) {

