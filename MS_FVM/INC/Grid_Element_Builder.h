#pragma once
#include "Grid_File_Type.h"
#include "Element.h"
#include "Profiler.h"
#include "Log.h"

#include <map>
#include <unordered_set>
#include <sstream>

template <size_t space_dimension>
struct Grid_Elements
{
	std::vector<Element<space_dimension>> cell_elements;
	std::vector<Element<space_dimension>> boundary_elements;
	std::vector<std::pair<Element<space_dimension>, Element<space_dimension>>> periodic_boundary_element_pairs;
	std::vector<Element<space_dimension>> inner_face_elements;
};


namespace ms {
	inline ElementType string_to_element_type(const std::string& str);
	template <typename T>
	std::vector<T> extract_by_index(const std::vector<T>& set, const std::vector<size_t>& indexes);
}


template <typename Grid_File_Type, size_t space_dimension>
class Grid_Element_Builder;


template <size_t space_dimension>
class Grid_Element_Builder<Gmsh, space_dimension>
{
	using Space_Vector_ = EuclideanVector<space_dimension>;

public:
	Grid_Element_Builder(void) = delete;

	static Grid_Elements<space_dimension> build_from_grid_file(const std::string& grid_file_name);

	//private: for test
	static Text read_about(std::ifstream& grid_file_stream, const std::string& target);
	static std::vector<Space_Vector_> make_node_datas(const Text& node_text);
	static Grid_Elements<space_dimension> make_elements(const Text& element_text, const Text& physical_name_text, const std::vector<Space_Vector_>& node_datas);
	static std::vector<Element<space_dimension>> make_inner_face_elements(const std::vector<Element<space_dimension>>& cell_elements, const std::vector<Element<space_dimension>>& boundary_elements, const std::vector<Element<space_dimension>>& periodic_boundary_elements);
	static std::vector<std::pair<Element<space_dimension>, Element<space_dimension>>> match_periodic_boundaries(std::vector<Element<space_dimension>>& periodic_boundary_elements);
};


//template definition part
template <size_t space_dimension>
Grid_Elements<space_dimension> Grid_Element_Builder<Gmsh, space_dimension>::build_from_grid_file(const std::string& grid_file_name) {
	SET_TIME_POINT;
	Log::content_ << "================================================================================\n";
	Log::content_ << "\t\t\t\t PreProcessing \n";
	Log::content_ << "================================================================================\n";	
	Log::print(); 
		

	SET_TIME_POINT;
	
	const auto grid_file_path = "RSC/Grid/" + grid_file_name + ".msh";

	std::ifstream grid_file_stream(grid_file_path);
	dynamic_require(grid_file_stream.is_open(), "fail to open " + grid_file_path);
	
	const auto node_text			= read_about(grid_file_stream, "Nodes");
	const auto node_datas			= make_node_datas(node_text);
	
	const auto element_text			= read_about(grid_file_stream, "Elements");	
	const auto physical_name_text	= read_about(grid_file_stream, "PhysicalNames");

	Log::content_ << std::left << std::setw(50) << "@ Read Grid File" << " ----------- " << GET_TIME_DURATION << "s\n\n";
	Log::print();

	return make_elements(element_text, physical_name_text, node_datas);
}

template <size_t space_dimension>
std::vector<EuclideanVector<space_dimension>> Grid_Element_Builder<Gmsh, space_dimension>::make_node_datas(const Text& node_text) {
	std::vector<Space_Vector_> node_datas;
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
Grid_Elements<space_dimension> Grid_Element_Builder<Gmsh, space_dimension>::make_elements(const Text& element_text, const Text& physical_name_text, const std::vector<Space_Vector_>& node_datas) {
	SET_TIME_POINT;

	std::map<index, ElementType> physical_group_index_to_element_type;
	for (const auto& physical_name_sentence : physical_name_text) {
		const char delimiter = ' ';
		const auto parsed_sentence_set = ms::parse(physical_name_sentence, delimiter);

		//const size_t dimension		= parsed_sentence_set[0].toValue<size_t>();
		const auto physical_group_index	= ms::string_to_value<index>(parsed_sentence_set[1]);
		const auto name					= ms::erase(parsed_sentence_set[2], "\"");
		const auto element_type			= ms::string_to_element_type(name);

		physical_group_index_to_element_type.emplace(physical_group_index, element_type);
	}

	std::vector<Element<space_dimension>> cell_elements;
	std::vector<Element<space_dimension>> boundary_elements;
	std::vector<Element<space_dimension>> periodic_boundary_elements;
	for (const auto& element_sentence : element_text) {
		const auto delimiter = ' ';
		const auto parsed_sentences = ms::parse(element_sentence, delimiter);

		auto value_set = ms::string_to_value_set<index>(parsed_sentences);

		//const auto index					= value_set[0];
		const auto figure_type_index		= value_set[1];
		//const auto tag_index				= value_set[2];
		const auto physical_gorup_index		= value_set[3];
		//const auto element_group_index	= value_set[4];

		//reference geometry
		const auto figure		= Gmsh::figure_type_index_to_element_figure(figure_type_index);
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

	const auto inner_face_elements = make_inner_face_elements(cell_elements, boundary_elements, periodic_boundary_elements);
	const auto periodic_boundary_element_pairs = match_periodic_boundaries(periodic_boundary_elements);

	Log::content_ << std::left << std::setw(50) << "@ Make Elements" << " ----------- " << GET_TIME_DURATION << "s\n";
	Log::content_ << "  " << std::setw(8) << cell_elements.size() << " cell \n";
	Log::content_ << "  " << std::setw(8) << boundary_elements.size() << " boundary\n";
	Log::content_ << "  " << std::setw(8) << periodic_boundary_element_pairs.size() << " periodic boundary pair\n";
	Log::content_ << "  " << std::setw(8) << inner_face_elements.size() << " inner face\n\n";
	Log::print();
	return { cell_elements, boundary_elements, periodic_boundary_element_pairs, inner_face_elements };
}

template<size_t space_dimension>
std::vector<Element<space_dimension>> Grid_Element_Builder<Gmsh, space_dimension>::make_inner_face_elements(const std::vector<Element<space_dimension>>& cell_elements, const std::vector<Element<space_dimension>>& boundary_elements, const std::vector<Element<space_dimension>>& periodic_boundary_elements) {
	//construct inner face elements
	std::map<std::vector<size_t>, Element<space_dimension>> vnode_indexes_to_inner_face_element;
	for (const auto& cell_element : cell_elements) {
		auto inner_face_elements = cell_element.make_inner_face_elements();
		for (auto& inner_face_element : inner_face_elements) {
			auto vnode_indexes = inner_face_element.vertex_node_indexes();
			std::sort(vnode_indexes.begin(), vnode_indexes.end());	//to ignore index order
			vnode_indexes_to_inner_face_element.emplace(std::move(vnode_indexes), std::move(inner_face_element));
		}
	}

	//erase constructed elements
	for (const auto& boundray_element : boundary_elements) {
		auto vnode_indexes = boundray_element.vertex_node_indexes();
		std::sort(vnode_indexes.begin(), vnode_indexes.end());	//to ignore index order
		const auto result = vnode_indexes_to_inner_face_element.erase(vnode_indexes);
		dynamic_require(result == 1, "boundary geometry should be one of inner face");
	}
	for (const auto& periodic_boundray_element : periodic_boundary_elements) {
		auto vnode_indexes = periodic_boundray_element.vertex_node_indexes();
		std::sort(vnode_indexes.begin(), vnode_indexes.end());	//to ignore index order
		const auto result = vnode_indexes_to_inner_face_element.erase(vnode_indexes);
		dynamic_require(result == 1, "periodic boundary geometry should be one of inner face");
	}

	// container change
	std::vector<Element<space_dimension>> inner_face_elements;

	const auto num_inner_face = vnode_indexes_to_inner_face_element.size();
	inner_face_elements.reserve(num_inner_face);

	for (auto&& [key, value] : vnode_indexes_to_inner_face_element)
		inner_face_elements.push_back(std::move(value));

	return inner_face_elements;
}

template<size_t space_dimension>
std::vector<std::pair<Element<space_dimension>, Element<space_dimension>>> Grid_Element_Builder<Gmsh, space_dimension>::match_periodic_boundaries(std::vector<Element<space_dimension>>& periodic_boundary_elements) {
	
	const auto num_periodic_element = periodic_boundary_elements.size();
	const auto num_pair = static_cast<size_t>(0.5 * num_periodic_element);

	std::unordered_set<size_t> matched_index;
	matched_index.reserve(num_periodic_element);

	std::vector<std::pair<Element<space_dimension>, Element<space_dimension>>> matched_periodic_element_pairs;
	matched_periodic_element_pairs.reserve(num_pair);

	for (size_t i = 0; i < num_periodic_element; ++i) {
		if (matched_index.find(i) != matched_index.end())
			continue;

		for (size_t j = i + 1; j < num_periodic_element; ++j) {
			if (matched_index.find(j) != matched_index.end())
				continue;

			auto& i_element = periodic_boundary_elements[i];
			auto& j_element = periodic_boundary_elements[j];

			if (i_element.is_periodic_pair(j_element)) {
				matched_periodic_element_pairs.push_back(std::make_pair(std::move(i_element), std::move(j_element)));
				matched_index.insert(i);
				matched_index.insert(j);
				break;
			}
		}
	}

	return matched_periodic_element_pairs;
}



template <size_t space_dimension>
Text Grid_Element_Builder<Gmsh, space_dimension>::read_about(std::ifstream& grid_file_stream, const std::string& target) {
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
