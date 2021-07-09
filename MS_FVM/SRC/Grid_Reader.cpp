#include "../INC/Grid_Reader.h"
namespace ms {
	//std::string to_string(const ElementType element_type) {
	//	switch (element_type) {
	//	case ElementType::cell:
	//		return "Cell";
	//	case ElementType::slip_wall_2D:
	//		return "SlipWall2D";
	//	case ElementType::supersonic_outlet_2D:
	//		return "SuperSonicOutlet2D";
	//	case ElementType::x_periodic:
	//		return "Xperiodic";
	//	case ElementType::y_periodic:
	//		return "Yperiodic";
	//	default:
	//		throw std::runtime_error("wrong element_type");
	//		return std::string();
	//	}
	//}

	ElementType string_to_element_type(const std::string& str) {
		if (ms::is_there_icase(str, "Unspecified"))
			return ElementType::cell;
		else if (ms::is_there_icase(str,"slip") && ms::is_there_icase(str, "wall") && ms::is_there_icase(str,"2D"))
			return ElementType::slip_wall_2D;
		else if (ms::is_there_icase(str, "SuperSonic") && ms::is_there_icase(str, "outlet") && ms::is_there_icase(str, "2D"))
			return ElementType::supersonic_outlet_2D;
		else if (ms::is_there_icase(str, "x")&& ms::is_there_icase(str, "periodic"))
			return ElementType::x_periodic;
		else if (ms::is_there_icase(str, "y") && ms::is_there_icase(str, "periodic"))
			return ElementType::y_periodic;
		else {
			throw std::runtime_error("wrong element_type");
			return ElementType::not_in_list;
		}
	}
}

Grid_Data Gmsh_Grid_Reader::read(const std::string& grid_file_path) {
	std::ifstream grid_file_stream(grid_file_path);
	dynamic_require(grid_file_stream.is_open(), "fail to open grid file!");

	const auto node_text			= Gmsh_Grid_Reader::read_about(grid_file_stream, "Nodes");
	const auto node_grid_data_set	= Gmsh_Grid_Reader::make_node_grid_data(node_text);

	const auto element_text			= Gmsh_Grid_Reader::read_about(grid_file_stream, "Elements");
	const auto physical_name_text	= Gmsh_Grid_Reader::read_about(grid_file_stream, "PhysicalNames");
	const auto element_data_set		= Gmsh_Grid_Reader::make_element_data(element_text, physical_name_text);

	return { node_grid_data_set, element_data_set[0], element_data_set[1], element_data_set[2] };
}

std::vector<Physical_Domain_Vector> Gmsh_Grid_Reader::make_node_grid_data(const Text& node_text) {
	std::vector<Physical_Domain_Vector> node_datas;
	node_datas.reserve(node_text.size());
	for (const auto& node_data : node_text) {
		const char delimiter = ' ';
		auto parsed_data_set = ms::parse(node_data, delimiter);

		//const auto node_index = ms::string_to_value<size_t>(parsed_data_set[0]);
		std::array<double, s_physical_domain_dimension> node_coords;
		for (size_t i = 0; i < s_physical_domain_dimension; ++i)
			node_coords[i] = ms::string_to_value<double>(parsed_data_set[i + 1]);

		node_datas.push_back(node_coords);
	}
	return node_datas;
}

std::array<std::vector<ElementGridData>, 3> Gmsh_Grid_Reader::make_element_data(const Text& element_text, const Text& physical_name_text){
	std::map<size_t, ElementType> physical_group_index_to_element_type;
	for (const auto& physical_name_sentence : physical_name_text) {
		const char delimiter = ' ';
		const auto parsed_sentence_set = ms::parse(physical_name_sentence, delimiter);

		//const size_t dimension = parsed_sentence_set[0].toValue<size_t>();
		const auto index = ms::string_to_value<size_t>(parsed_sentence_set[1]);
		const auto name = ms::erase(parsed_sentence_set[2], "\"");
		const auto element_type = ms::string_to_element_type(name);

		physical_group_index_to_element_type.emplace(index, element_type);
	}

	std::vector<ElementGridData> cell_data;
	std::vector<ElementGridData> boundary_face_data;
	std::vector<ElementGridData> periodic_face_data;
	for (const auto& element_sentence : element_text) {
		const auto delimiter = ' ';
		auto parsed_sentence_set = ms::parse(element_sentence, delimiter);

		auto value_set = ms::string_to_value_set<size_t>(parsed_sentence_set);

		const auto index					= value_set[0];
		const auto figure_type_index		= value_set[1];
		//const auto num_tag				= value_set[2];
		const auto physical_gorup_index		= value_set[3];
		//const auto element_group_index	= value_set[4];

		constexpr size_t num_index = 5;
		value_set.erase(value_set.begin(), value_set.begin() + num_index);

		const auto figure				= Gmsh_Grid_Reader::figure_type_index_to_element_figure(figure_type_index);
		const auto figure_order			= Gmsh_Grid_Reader::figure_type_index_to_figure_order(figure_type_index);
		const auto type					= physical_group_index_to_element_type.at(physical_gorup_index);
		auto& consisting_node_index_set = value_set;

		ElementGridData element_data = { index, figure, figure_order, type, std::move(consisting_node_index_set) };

		switch (type) {
		case ElementType::cell:
			cell_data.emplace_back(std::move(element_data));
			break;
		case ElementType::x_periodic:
		case ElementType::y_periodic:
			periodic_face_data.emplace_back(std::move(element_data));
			break;
		default:
			boundary_face_data.emplace_back(std::move(element_data));
			break;
		}
	}

	return { cell_data, boundary_face_data, periodic_face_data };
}

size_t Gmsh_Grid_Reader::figure_type_index_to_figure_order(const size_t element_type_indx) {
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

ElementFigure Gmsh_Grid_Reader::figure_type_index_to_element_figure(const size_t element_type_index) {
	switch (static_cast<GmshFigureType>(element_type_index)) {
	case GmshFigureType::POINT:			return ElementFigure::point;
	case GmshFigureType::LINE_P1:
	case GmshFigureType::LINE_P2:
	case GmshFigureType::LINE_P3:
	case GmshFigureType::LINE_P4:
	case GmshFigureType::LINE_P5:
	case GmshFigureType::LINE_P6:		return ElementFigure::line;
	case GmshFigureType::TRIS_P1:
	case GmshFigureType::TRIS_P2:
	case GmshFigureType::TRIS_P3:
	case GmshFigureType::TRIS_P4:
	case GmshFigureType::TRIS_P5:		return ElementFigure::triangle;
	case GmshFigureType::QUAD_P1:
	case GmshFigureType::QUAD_P2:
	case GmshFigureType::QUAD_P3:
	case GmshFigureType::QUAD_P4:
	case GmshFigureType::QUAD_P5:
	case GmshFigureType::QUAD_P6:		return ElementFigure::quadrilateral;
	case GmshFigureType::TETS_P1:
	case GmshFigureType::TETS_P2:
	case GmshFigureType::TETS_P3:
	case GmshFigureType::TETS_P4:
	case GmshFigureType::TETS_P5:		return ElementFigure::tetrahedral;
	case GmshFigureType::HEXA_P1:
	case GmshFigureType::HEXA_P2:
	case GmshFigureType::HEXA_P3:
	case GmshFigureType::HEXA_P4:
	case GmshFigureType::HEXA_P5:		return ElementFigure::hexahedral;
	case GmshFigureType::PRIS_P1:
	case GmshFigureType::PRIS_P2:
	case GmshFigureType::PRIS_P3:
	case GmshFigureType::PRIS_P4:
	case GmshFigureType::PRIS_P5:		return ElementFigure::prism;
	case GmshFigureType::PYRA_P1:
	case GmshFigureType::PYRA_P2:
	case GmshFigureType::PYRA_P3:
	case GmshFigureType::PYRA_P4:		return ElementFigure::pyramid;
	default:
		throw std::runtime_error("invalid element type index");
		return ElementFigure::not_in_list;
	}
}

Text Gmsh_Grid_Reader::read_about(std::ifstream& grid_file_stream, const std::string& target) {
	const auto target_str = "$" + target;

	std::string tmp_str;
	while (std::getline(grid_file_stream, tmp_str)) {
		if (tmp_str.find(target_str) != std::string::npos) {
			std::getline(grid_file_stream, tmp_str);
			break;
		}
	}

	const auto num_data = ms::string_to_value<size_t>(tmp_str);
	return Text(grid_file_stream, num_data);
}