#include "../INC/Grid_File_Type.h"

order Gmsh::figure_type_index_to_figure_order(const index element_type_indx) {
	switch (static_cast<GmshFigureType>(element_type_indx)) {
	case GmshFigureType::POINT:			return 0;
	case GmshFigureType::LINE_P1:
	case GmshFigureType::TRIS_P1:
	case GmshFigureType::QUAD_P1:
	case GmshFigureType::TETS_P1:
	case GmshFigureType::HEXA_P1:
	case GmshFigureType::PRIS_P1:
	case GmshFigureType::PYRA_P1:		return 1;
	case GmshFigureType::LINE_P2:
	case GmshFigureType::TRIS_P2:
	case GmshFigureType::QUAD_P2:
	case GmshFigureType::TETS_P2:
	case GmshFigureType::HEXA_P2:
	case GmshFigureType::PRIS_P2:
	case GmshFigureType::PYRA_P2:		return 2;
	case GmshFigureType::LINE_P3:
	case GmshFigureType::TRIS_P3:
	case GmshFigureType::QUAD_P3:
	case GmshFigureType::TETS_P3:
	case GmshFigureType::HEXA_P3:
	case GmshFigureType::PRIS_P3:
	case GmshFigureType::PYRA_P3:		return 3;
	case GmshFigureType::LINE_P4:
	case GmshFigureType::TRIS_P4:
	case GmshFigureType::QUAD_P4:
	case GmshFigureType::TETS_P4:
	case GmshFigureType::HEXA_P4:
	case GmshFigureType::PRIS_P4:
	case GmshFigureType::PYRA_P4:		return 4;
	case GmshFigureType::LINE_P5:
	case GmshFigureType::TRIS_P5:
	case GmshFigureType::QUAD_P5:
	case GmshFigureType::TETS_P5:
	case GmshFigureType::HEXA_P5:
	case GmshFigureType::PRIS_P5:		return 5;
	case GmshFigureType::LINE_P6:
	case GmshFigureType::QUAD_P6:		return 6;
	default:
		throw std::runtime_error("invalid element type index");
		return NULL;
	}
}


Figure Gmsh::figure_type_index_to_element_figure(const index element_type_index) {
	switch (static_cast<GmshFigureType>(element_type_index)) {
	case GmshFigureType::POINT:			return Figure::point;
	case GmshFigureType::LINE_P1:
	case GmshFigureType::LINE_P2:
	case GmshFigureType::LINE_P3:
	case GmshFigureType::LINE_P4:
	case GmshFigureType::LINE_P5:
	case GmshFigureType::LINE_P6:		return Figure::line;
	case GmshFigureType::TRIS_P1:
	case GmshFigureType::TRIS_P2:
	case GmshFigureType::TRIS_P3:
	case GmshFigureType::TRIS_P4:
	case GmshFigureType::TRIS_P5:		return Figure::triangle;
	case GmshFigureType::QUAD_P1:
	case GmshFigureType::QUAD_P2:
	case GmshFigureType::QUAD_P3:
	case GmshFigureType::QUAD_P4:
	case GmshFigureType::QUAD_P5:
	case GmshFigureType::QUAD_P6:		return Figure::quadrilateral;
	case GmshFigureType::TETS_P1:
	case GmshFigureType::TETS_P2:
	case GmshFigureType::TETS_P3:
	case GmshFigureType::TETS_P4:
	case GmshFigureType::TETS_P5:		return Figure::tetrahedral;
	case GmshFigureType::HEXA_P1:
	case GmshFigureType::HEXA_P2:
	case GmshFigureType::HEXA_P3:
	case GmshFigureType::HEXA_P4:
	case GmshFigureType::HEXA_P5:		return Figure::hexahedral;
	case GmshFigureType::PRIS_P1:
	case GmshFigureType::PRIS_P2:
	case GmshFigureType::PRIS_P3:
	case GmshFigureType::PRIS_P4:
	case GmshFigureType::PRIS_P5:		return Figure::prism;
	case GmshFigureType::PYRA_P1:
	case GmshFigureType::PYRA_P2:
	case GmshFigureType::PYRA_P3:
	case GmshFigureType::PYRA_P4:		return Figure::pyramid;
	default:
		throw std::runtime_error("invalid element type index");
		return Figure::not_in_list;
	}
}