#pragma once
#include "Grid_File_Type.h"
#include "Text.h"
#include "Element.h"
#include "Profiler.h"

#include <map>


template <size_t space_dimension>
struct Grid_Raw_Data
{
	std::vector<Element<space_dimension>> cell_elements;
	std::vector<Element<space_dimension>> boundary_elements;
	std::vector<Element<space_dimension>> periodic_boundary_elements;
};


namespace ms {
	inline ElementType string_to_element_type(const std::string& str);
}


template <typename Grid_File_Type, size_t space_dimension>
class Grid_File_Convertor;


template <size_t space_dimension>
class Grid_File_Convertor<Gmsh, space_dimension>
{
	using Space_Vector = EuclideanVector<space_dimension>;

public:
	Grid_File_Convertor(void) = delete;

	static Grid_Raw_Data<space_dimension> convert(const std::string& grid_file_name);

	//private: for test
	static Text read_about(std::ifstream& grid_file_stream, const std::string& target);
	static std::vector<Space_Vector> make_node_datas(const Text& node_text);
	static std::array<std::vector<Element<space_dimension>>, 3> make_element_datas(const Text& element_text, const Text& physical_name_text, const std::vector<Space_Vector>& node_datas);
};


//template definition part
template <size_t space_dimension>
Grid_Raw_Data<space_dimension> Grid_File_Convertor<Gmsh, space_dimension>::convert(const std::string& grid_file_name) {
	std::cout << std::left;
	std::cout << "============================================================\n";
	std::cout << "\t Grid construction start!\n";
	std::cout << "============================================================\n";
	SET_TIME_POINT;

	const auto grid_file_path = "RSC/Grid/" + grid_file_name + ".msh";

	std::ifstream grid_file_stream(grid_file_path);
	dynamic_require(grid_file_stream.is_open(), "fail to open grid file!");
	
	const auto node_text			= read_about(grid_file_stream, "Nodes");
	const auto node_datas			= make_node_datas(node_text);
	
	const auto element_text			= read_about(grid_file_stream, "Elements");
	const auto physical_name_text	= read_about(grid_file_stream, "PhysicalNames");
	const auto element_datas		= make_element_datas(element_text, physical_name_text, node_datas);

	return { element_datas[0], element_datas[1], element_datas[2] };
}

template <size_t space_dimension>
std::vector<EuclideanVector<space_dimension>> Grid_File_Convertor<Gmsh, space_dimension>::make_node_datas(const Text& node_text) {
	std::vector<Space_Vector> node_datas;
	node_datas.reserve(node_text.size());
	for (const auto& node_data : node_text) {
		const char delimiter = ' ';
		auto parsed_data_set = ms::parse(node_data, delimiter);

		//const auto node_index = ms::string_to_value<size_t>(parsed_data_set[0]);
		std::array<double, space_dimension> node_coords;
		for (size_t i = 0; i < space_dimension; ++i)
			node_coords[i] = ms::string_to_value<double>(parsed_data_set[i + 1]);

		node_datas.push_back(node_coords);
	}
	return node_datas;
}

template <size_t space_dimension>
std::array<std::vector<Element<space_dimension>>, 3> Grid_File_Convertor<Gmsh, space_dimension>::make_element_datas(const Text& element_text, const Text& physical_name_text, const std::vector<Space_Vector>& node_datas) {
	SET_TIME_POINT;

	std::map<size_t, ElementType> physical_group_index_to_element_type;
	for (const auto& physical_name_sentence : physical_name_text) {
		const char delimiter = ' ';
		const auto parsed_sentence_set = ms::parse(physical_name_sentence, delimiter);

		//const size_t dimension		= parsed_sentence_set[0].toValue<size_t>();
		const auto index				= ms::string_to_value<size_t>(parsed_sentence_set[1]);
		const auto name					= ms::erase(parsed_sentence_set[2], "\"");
		const auto element_type			= ms::string_to_element_type(name);

		physical_group_index_to_element_type.emplace(index, element_type);
	}

	std::vector<Element<space_dimension>> cell_elements;
	std::vector<Element<space_dimension>> boundary_elements;
	std::vector<Element<space_dimension>> periodic_boundary_elements;
	for (const auto& element_sentence : element_text) {
		const auto delimiter = ' ';
		const auto parsed_sentences = ms::parse(element_sentence, delimiter);

		auto value_set = ms::string_to_value_set<size_t>(parsed_sentences);

		//const auto index					= value_set[0];
		const auto figure_type_index		= value_set[1];
		//const auto tag_index				= value_set[2];
		const auto physical_gorup_index		= value_set[3];
		//const auto element_group_index	= value_set[4];

		//reference geometry
		const auto figure = Gmsh::figure_type_index_to_element_figure(figure_type_index);
		const auto figure_order = Gmsh::figure_type_index_to_figure_order(figure_type_index);
		auto reference_geometry = ReferenceGeometry(figure, figure_order);

		//geometry
		constexpr size_t num_index = 5;
		value_set.erase(value_set.begin(), value_set.begin() + num_index);

		const auto num_nodes = value_set.size();
		std::vector<size_t> node_indexes(num_nodes);
		for (size_t i = 0; i < num_nodes; ++i)
			node_indexes[i] = value_set[i] - 1;		//Gmsh node index start with 1

		auto nodes = ms::extract_by_index(node_datas, node_indexes);
		Geometry geometry(reference_geometry, std::move(nodes));

		//element
		const auto type	= physical_group_index_to_element_type.at(physical_gorup_index);
		Element element(type, std::move(geometry), std::move(node_indexes));
		

		switch (type) {
		case ElementType::cell:
			cell_elements.emplace_back(std::move(element));
			break;
		case ElementType::x_periodic:
		case ElementType::y_periodic:
			periodic_boundary_elements.emplace_back(std::move(element));
			break;
		default:
			boundary_elements.emplace_back(std::move(element));
			break;
		}
	}

	std::cout << "make elements \t\t\t ellapsed " << std::setw(10) << GET_TIME_DURATION << "s\n";
	std::cout << "num cell :\t\t\t" << cell_elements.size() << "\n";
	std::cout << "num boundary : \t\t" << boundary_elements.size() << "\n";
	std::cout << "num periodic boundary : \t" << periodic_boundary_elements.size() << "\n";

	return { cell_elements, boundary_elements, periodic_boundary_elements };
}




template <size_t space_dimension>
Text Grid_File_Convertor<Gmsh, space_dimension>::read_about(std::ifstream& grid_file_stream, const std::string& target) {
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


//inline function definition
namespace ms {
	inline ElementType string_to_element_type(const std::string& str) {
		if (ms::is_there_icase(str, "Unspecified"))
			return ElementType::cell;
		else if (ms::is_there_icase(str, "slip") && ms::is_there_icase(str, "wall") && ms::is_there_icase(str, "2D"))
			return ElementType::slip_wall_2D;
		else if (ms::is_there_icase(str, "SuperSonic") && ms::is_there_icase(str, "outlet") && ms::is_there_icase(str, "2D"))
			return ElementType::supersonic_outlet_2D;
		else if (ms::is_there_icase(str, "x") && ms::is_there_icase(str, "periodic"))
			return ElementType::x_periodic;
		else if (ms::is_there_icase(str, "y") && ms::is_there_icase(str, "periodic"))
			return ElementType::y_periodic;
		else {
			throw std::runtime_error("wrong element_type");
			return ElementType::not_in_list;
		}
	}

	template <typename T>
	std::vector<T> extract_by_index(const std::vector<T>& set, const std::vector<size_t>& indexes) {
		const auto num_extracted_value = indexes.size();
		std::vector<T> extracted_values(num_extracted_value);
		for (size_t i = 0; i < num_extracted_value; ++i)
			extracted_values[i] = set[indexes[i]]; // index start with 1

		return extracted_values;
	}
}
