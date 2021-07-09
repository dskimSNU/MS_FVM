#pragma once
#include "EuclideanVector.h"
#include "Text.h"
#include "../RSC/Setting.h"
#include <map>


using Physical_Domain_Vector = EuclideanVector<PHYSICAL_DOMAIN_DIMENSION>;


enum class ElementFigure
{	
	point,	
	line,
	triangle,		quadrilateral,
	tetrahedral,	hexahedral,		prism,	pyramid,
	not_in_list
};


enum class ElementType
{
	cell,
	slip_wall_2D,
	supersonic_outlet_2D,
	x_periodic, y_periodic,
	not_in_list
};

struct ElementGridData
{
	size_t index;
	ElementFigure figure;
	size_t figure_order;
	ElementType type;
	std::vector<size_t> node_indexes;
};


struct Grid_Data
{
	std::vector<Physical_Domain_Vector> node_grid_datas;
	std::vector<ElementGridData> cell_grid_datas;
	std::vector<ElementGridData> boundary_grid_datas;
	std::vector<ElementGridData> periodic_boundary_grid_datas;
};


namespace ms {
	//std::string to_string(const ElementType element_type);
	ElementType string_to_element_type(const std::string& str);
}


enum class GmshFigureType
{
	POINT = 0,
	LINE_P1 = 1, LINE_P2 = 8,  LINE_P3 = 26,  LINE_P4 = 27, LINE_P5 = 28, LINE_P6 = 62,
	TRIS_P1 = 2, TRIS_P2 = 9,  TRIS_P3 = 21,  TRIS_P4 = 23, TRIS_P5 = 25,
	QUAD_P1 = 3, QUAD_P2 = 10, QUAD_P3 = 36,  QUAD_P4 = 37, QUAD_P5 = 38, QUAD_P6 = 47,
	TETS_P1 = 4, TETS_P2 = 11, TETS_P3 = 29,  TETS_P4 = 30, TETS_P5 = 31,
	HEXA_P1 = 5, HEXA_P2 = 12, HEXA_P3 = 92,  HEXA_P4 = 93, HEXA_P5 = 94,
	PRIS_P1 = 6, PRIS_P2 = 13, PRIS_P3 = 90,  PRIS_P4 = 91, PRIS_P5 = 106,
	PYRA_P1 = 7, PYRA_P2 = 14, PYRA_P3 = 118, PYRA_P4 = 119
};


class Gmsh_Grid_Reader
{
public:
	Gmsh_Grid_Reader(void) = delete;
	static Grid_Data read(const std::string& grid_file_name);

//private: for test
	static Text read_about(std::ifstream& grid_file_stream, const std::string& target);
	static std::vector<Physical_Domain_Vector> make_node_grid_data(const Text& node_text);
	static std::array<std::vector<ElementGridData>, 3> make_element_data(const Text& element_text, const Text& physical_name_text);
	static size_t figure_type_index_to_figure_order(const size_t figure_type_indx);
	static ElementFigure figure_type_index_to_element_figure(const size_t figure_type_index);
};