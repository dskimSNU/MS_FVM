#pragma once
#include "Grid_File_Type.h"
#include "Text.h"

#include <map>


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
	Figure figure;
	size_t figure_order;
	ElementType type;
	std::vector<size_t> node_indexes;
};


template <size_t dim>
struct Grid_Data
{
	std::vector<EuclideanVector<dim>> node_datas;
	std::vector<ElementGridData> cell_datas;
	std::vector<ElementGridData> boundary_datas;
	std::vector<ElementGridData> periodic_boundary_datas;
};


namespace ms {
	inline ElementType string_to_element_type(const std::string& str);
}


template <typename Grid_File_Type, size_t dim>
class Grid_File_To_Data;


template <size_t dim>
class Grid_File_To_Data<Gmsh, dim>
{
	using Space_Vector = EuclideanVector<dim>;

public:
	Grid_File_To_Data(void) = delete;

	static Grid_Data<dim> convert(const std::string& grid_file_name);

	//private: for test
	static Text read_about(std::ifstream& grid_file_stream, const std::string& target);
	static std::vector<Space_Vector> make_node_grid_data(const Text& node_text);
	static std::array<std::vector<ElementGridData>, 3> make_element_data(const Text& element_text, const Text& physical_name_text);

};


//template definition part
template <size_t dim>
Grid_Data<dim> Grid_File_To_Data<Gmsh, dim>::convert(const std::string& grid_file_path) {
	std::ifstream grid_file_stream(grid_file_path);
	dynamic_require(grid_file_stream.is_open(), "fail to open grid file!");
	
	const auto node_text			= read_about(grid_file_stream, "Nodes");
	const auto node_grid_datas		= make_node_grid_data(node_text);
	
	const auto element_text			= read_about(grid_file_stream, "Elements");
	const auto physical_name_text	= read_about(grid_file_stream, "PhysicalNames");
	const auto element_datas		= make_element_data(element_text, physical_name_text);

	return { node_grid_datas, element_datas[0], element_datas[1], element_datas[2] };
}

template <size_t dim>
std::vector<EuclideanVector<dim>> Grid_File_To_Data<Gmsh, dim>::make_node_grid_data(const Text& node_text) {

	std::vector<Space_Vector> node_datas;
	node_datas.reserve(node_text.size());
	for (const auto& node_data : node_text) {
		const char delimiter = ' ';
		auto parsed_data_set = ms::parse(node_data, delimiter);

		//const auto node_index = ms::string_to_value<size_t>(parsed_data_set[0]);
		std::array<double, dim> node_coords;
		for (size_t i = 0; i < dim; ++i)
			node_coords[i] = ms::string_to_value<double>(parsed_data_set[i + 1]);

		node_datas.push_back(node_coords);
	}
	return node_datas;
}

template <size_t dim>
std::array<std::vector<ElementGridData>, 3> Grid_File_To_Data<Gmsh, dim>::make_element_data(const Text& element_text, const Text& physical_name_text) {
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

	std::vector<ElementGridData> cell_data;
	std::vector<ElementGridData> boundary_face_data;
	std::vector<ElementGridData> periodic_face_data;
	for (const auto& element_sentence : element_text) {
		const auto delimiter = ' ';
		const auto parsed_sentences = ms::parse(element_sentence, delimiter);

		auto value_set = ms::string_to_value_set<size_t>(parsed_sentences);

		const auto index					= value_set[0];
		const auto figure_type_index		= value_set[1];
		//const auto tag_index				= value_set[2];
		const auto physical_gorup_index		= value_set[3];
		//const auto element_group_index	= value_set[4];

		constexpr size_t num_index = 5;
		value_set.erase(value_set.begin(), value_set.begin() + num_index);
		const auto figure				= Gmsh::figure_type_index_to_element_figure(figure_type_index);
		const auto figure_order			= Gmsh::figure_type_index_to_figure_order(figure_type_index);
		const auto type					= physical_group_index_to_element_type.at(physical_gorup_index);
		auto& consisting_node_indexes	= value_set;

		ElementGridData element_data = { index, figure, figure_order, type, std::move(consisting_node_indexes) };

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



template <size_t dim>
Text Grid_File_To_Data<Gmsh, dim>::read_about(std::ifstream& grid_file_stream, const std::string& target) {
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


//inline free function definition
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
}
