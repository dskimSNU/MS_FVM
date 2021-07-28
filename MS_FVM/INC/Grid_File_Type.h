#pragma once
#include "Element.h"

#include <type_traits>

using index = unsigned int;
using order = unsigned short;


class GFT {}; //Grid File Type


enum class GmshFigureType
{
	POINT = 0,
	LINE_P1 = 1, LINE_P2 = 8, LINE_P3 = 26, LINE_P4 = 27, LINE_P5 = 28, LINE_P6 = 62,
	TRIS_P1 = 2, TRIS_P2 = 9, TRIS_P3 = 21, TRIS_P4 = 23, TRIS_P5 = 25,
	QUAD_P1 = 3, QUAD_P2 = 10, QUAD_P3 = 36, QUAD_P4 = 37, QUAD_P5 = 38, QUAD_P6 = 47,
	TETS_P1 = 4, TETS_P2 = 11, TETS_P3 = 29, TETS_P4 = 30, TETS_P5 = 31,
	HEXA_P1 = 5, HEXA_P2 = 12, HEXA_P3 = 92, HEXA_P4 = 93, HEXA_P5 = 94,
	PRIS_P1 = 6, PRIS_P2 = 13, PRIS_P3 = 90, PRIS_P4 = 91, PRIS_P5 = 106,
	PYRA_P1 = 7, PYRA_P2 = 14, PYRA_P3 = 118, PYRA_P4 = 119
};


class Gmsh : public GFT {
public:
	static order figure_type_index_to_figure_order(const index figure_type_index);
	static Figure figure_type_index_to_element_figure(const index figure_type_index);
};


namespace ms {
	template <typename T>
	inline constexpr bool is_grid_file_type = std::is_base_of_v<GFT, T>;
}