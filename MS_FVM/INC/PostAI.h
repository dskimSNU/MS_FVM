#pragma once
#include "Grid_Builder.h"
#include "Governing_Equation.h"
#include "Inital_Condition.h"
#include "Setting.h"

class PostAI
{
//private: //for test
public:
	inline static std::string path_;	//inline static vs static?
	inline static size_t num_data_;
	inline static Text comment;			//comment 추가
	//inline static std::vector<std::vector<size_t>> set_of_index_orders_; //김동석이 숙제


	inline static std::vector<Text> ai_data_text_set_;
	inline static std::vector<std::vector<size_t>> vertex_share_cell_indexes_set_;
	inline static std::vector<size_t> target_cell_indexes_;

private:
	PostAI(void) = delete;

public:
	static void set_path(const std::string& path);

	template <size_t space_dimension>
	static void intialize(const Grid<space_dimension>& grid);

	template <size_t num_equation>
	static void record_solution_datas(const std::vector<EuclideanVector<num_equation>>& solutions, const std::vector<Dynamic_Matrix_>& solution_gradients);	

	template <size_t num_equation>
	static void record_limiting_value(const size_t index, const std::array<double, num_equation>& limiting_value);


	static void post(void);

//private: //for test
	template <size_t space_dimension>
	static auto calculate_face_share_cell_indexes_set(const Grid<space_dimension>& grid);

	template <size_t space_dimension>
	static auto calculate_vertex_nodes_coordinate_string_set(const Grid<space_dimension>& grid);

	template <size_t num_equation>
	static auto convert_to_solution_strings(const std::vector<EuclideanVector<num_equation>>& solutions);	

	static std::vector<std::string> convert_to_solution_gradient_strings(const std::vector<Dynamic_Matrix_>& solution_gradients);
};


//template definition part
template <size_t space_dimension>
void PostAI::intialize(const Grid<space_dimension>& grid) {
#ifdef POST_AI_DATA

	const auto& vnode_index_to_share_cell_indexes = grid.connectivity.vnode_index_to_share_cell_indexes;
	const auto& cell_elements = grid.elements.cell_elements;
	num_data_ = cell_elements.size();

	vertex_share_cell_indexes_set_.reserve(num_data_);
	ai_data_text_set_.resize(num_data_);
	target_cell_indexes_.reserve(num_data_);
	//set_of_index_orders_.reserve(num_data_);

	auto file_path = path_ + "AI_Solver_Data_1.txt";	//추가
	const auto parsed_path = ms::parse(path_, '/');		
	comment << "***********************************\n"
		"***********************************\n"
		"Grid : " + parsed_path[7] + "\n" +
		"GE : " + parsed_path[4] + "\n" +
		"IC : " + parsed_path[5] + "\n"
		"***********************************\n"
		"***********************************\n\n";
	comment.add_write(file_path);


	const auto face_share_cell_indexes_set = calculate_face_share_cell_indexes_set(grid);
	const auto vnodes_coordinate_string_set = calculate_vertex_nodes_coordinate_string_set(grid);

	for (size_t i = 0; i < num_data_; ++i) {
		const auto& cell_element = cell_elements[i];
		const auto& cell_geometry = cell_element.geometry_;
		

		//vertex share cell indexes temp
		std::set<size_t> vertex_share_cell_indexes_temp;

		const auto vnode_indexes = cell_element.vertex_node_indexes();
		for (const auto vnode_index : vnode_indexes) {
			const auto& share_cell_indexes = vnode_index_to_share_cell_indexes.at(vnode_index);
			vertex_share_cell_indexes_temp.insert(share_cell_indexes.begin(), share_cell_indexes.end());
		}

		////make face share cell indexes
		//std::vector<size_t> face_share_cell_indexes; //동석이 숙제

		//const auto face_vnode_indexes_set = cell_element.face_vertex_node_indexes_set();
		//for (const auto& face_vnode_indexes : face_vnode_indexes_set) {
		//	//only work at 2D
		//	std::vector<size_t> temp;

		//	dynamic_require(face_vnode_indexes.size() == 2, "This algorithm only work at 2D");
		//	const auto vnode_index0 = face_vnode_indexes[0];
		//	const auto vnode_index1 = face_vnode_indexes[1];

		//	const auto& set0 = vnode_index_to_share_cell_indexes.at(vnode_index0);
		//	const auto& set1 = vnode_index_to_share_cell_indexes.at(vnode_index1);

		//	std::set_intersection(set0.begin(), set0.end(), set1.begin(), set1.end(), std::back_inserter(temp));
		//	dynamic_require(temp.size() == 2, "face share cell should be unique");

		//	const auto& target_cell_index = i;
		//	const auto target_cell_index_iter = std::find(temp.begin(), temp.end(), target_cell_index);			
		//	temp.erase(target_cell_index_iter);
		//	
		//	face_share_cell_indexes.push_back(temp.front());
		//}

		////make cross vnode cell indexes
		//std::set<size_t> set_of_cell_face_share_cell_indexes(face_share_cell_indexes.begin(), face_share_cell_indexes.end());

		//std::vector<size_t> cross_vnode_cell_indexes;
		//
		//dynamic_require(vnode_indexes.size() == 4, "This algorithm only work for quadrilateral element");
		//for (const auto vnode_index : vnode_indexes) {
		//	const auto& share_cell_indexes = vnode_index_to_share_cell_indexes.at(vnode_index);

		//	std::vector<size_t> temp; // 동석이 숙제
		//	std::set_difference(share_cell_indexes.begin(), share_cell_indexes.end(), set_of_cell_face_share_cell_indexes.begin(), set_of_cell_face_share_cell_indexes.end(), std::back_inserter(temp));
		//	dynamic_require(temp.size() == 2, "This algorithm only work for quadrilateral element");

		//	const auto& target_cell_index = i;
		//	const auto target_cell_index_iter = std::find(temp.begin(), temp.end(), target_cell_index);
		//	temp.erase(target_cell_index_iter);

		//	cross_vnode_cell_indexes.push_back(temp.front());
		//}

		////make index orders
		//auto& index_orders = set_of_index_orders_[i];
		//for (size_t j = 0; j < 4; ++j) {
		//	if (j == 3) {
		//		index_orders.push_back(face_share_cell_indexes[j]);
		//		index_orders.push_back(cross_vnode_cell_indexes.front());
		//	}

		//	index_orders.push_back(face_share_cell_indexes[j]);
		//	index_orders.push_back(cross_vnode_cell_indexes[j + 1]);
		//}



		//chunk edge connectivities //quad3에서는 제대로 작동하지 않는 algorithm
		std::set<std::set<size_t>> face_share_cell_index_pairs;

		for (const auto chunk_cell_index : vertex_share_cell_indexes_temp) {
			const auto& face_share_cell_indexes = face_share_cell_indexes_set.at(chunk_cell_index);

			std::vector<size_t> face_share_cell_indexes_in_chunk;
			std::set_intersection(vertex_share_cell_indexes_temp.begin(), vertex_share_cell_indexes_temp.end(), face_share_cell_indexes.begin(), face_share_cell_indexes.end(), std::back_inserter(face_share_cell_indexes_in_chunk));

			for (const auto face_share_cell_index_in_chunk : face_share_cell_indexes_in_chunk)
				face_share_cell_index_pairs.insert({ chunk_cell_index, face_share_cell_index_in_chunk });
		}

		//vertex_share_cell_indexes
		vertex_share_cell_indexes_temp.erase(i);

		std::vector<size_t> vertex_share_cell_indexes;
		vertex_share_cell_indexes.push_back(i);
		vertex_share_cell_indexes.insert(vertex_share_cell_indexes.end(), vertex_share_cell_indexes_temp.begin(), vertex_share_cell_indexes_temp.end());

		// header string
		ai_data_text_set_[i] << "temporary header sentence"; // 나중에 ## 들어올 자리

		// node number string
		const auto num_node = vertex_share_cell_indexes.size();
		ai_data_text_set_[i] << "@nodeNumber\n" + std::to_string(num_node);	

		// node index order string
		std::string node_index_order_string = "@nodeIndexOrder\n";
		for (const auto& node_index : vertex_share_cell_indexes)
			node_index_order_string += std::to_string(node_index) + "\t";
		ai_data_text_set_[i] << std::move(node_index_order_string);	//move : 형변환? 다시 공부

		//edge number string
		const auto num_edge = face_share_cell_index_pairs.size();
		ai_data_text_set_[i] << "@edgeNumber\n" + std::to_string(num_edge);

		//connectivity string
		std::string node_connectivity_string = "@connectivity\n";
		for (const auto& chunk_edge_connectivity : face_share_cell_index_pairs) {
			for (const auto& node_index : chunk_edge_connectivity)
				node_connectivity_string += std::to_string(node_index) + "\t";
			node_connectivity_string += "\n";
		}
		node_connectivity_string.pop_back();
		ai_data_text_set_[i] << std::move(node_connectivity_string);

		//cell coords string
		std::string cell_coords_string = "@cellCoords\n";
		for (const auto vertex_share_cell_index : vertex_share_cell_indexes)
			cell_coords_string += vnodes_coordinate_string_set[vertex_share_cell_index];
		cell_coords_string.pop_back();
		ai_data_text_set_[i] << std::move(cell_coords_string);

		//vertex_share_cell_indexes
		vertex_share_cell_indexes_set_.push_back(std::move(vertex_share_cell_indexes));
	}

#endif
}


template <size_t num_equation>
void PostAI::record_solution_datas(const std::vector<EuclideanVector<num_equation>>& solutions, const std::vector<Dynamic_Matrix_>& solution_gradients) {
#ifdef POST_AI_DATA

	dynamic_require(num_data_ == solutions.size(),			"number of solution should be same with number of data");
	dynamic_require(num_data_ == solution_gradients.size(), "number of solution gradient should be same with number of data");

	const auto solution_strings = convert_to_solution_strings(solutions);
	const auto solution_gradient_strings = convert_to_solution_gradient_strings(solution_gradients);

	std::string cell_average_string;
	std::string cell_gradient_string;
	for (size_t i = 0; i < num_data_; ++i) {
		const auto& vertex_share_cell_indexes = vertex_share_cell_indexes_set_[i];
		const auto num_vertex_share_cell = vertex_share_cell_indexes.size();


		std::vector<double> vertex_share_cell_solutions(num_vertex_share_cell);
		for (size_t i = 0; i < num_vertex_share_cell; ++i)
			vertex_share_cell_solutions[i] = solutions[vertex_share_cell_indexes[i]][0];

		const auto min_solution = *std::min_element(vertex_share_cell_solutions.begin(), vertex_share_cell_solutions.end());
		const auto max_solution = *std::max_element(vertex_share_cell_solutions.begin(), vertex_share_cell_solutions.end());
		const auto solution_diff = max_solution - min_solution;

		if (solution_diff < 0.01)
			continue;
			
		target_cell_indexes_.push_back(i);

		cell_average_string = "@cellAverage\n";
		cell_gradient_string = "@cellGradient\n";
		for (const auto vertex_share_cell_index : vertex_share_cell_indexes) {
			cell_average_string += solution_strings[vertex_share_cell_index] + "\n";
			cell_gradient_string += solution_gradient_strings[vertex_share_cell_index] + "\n";
		}
		cell_average_string.pop_back();
		cell_gradient_string.pop_back();

		ai_data_text_set_[i] << std::move(cell_average_string) << std::move(cell_gradient_string);
	}

#endif
}

template <size_t num_equation>
void PostAI::record_limiting_value(const size_t index, const std::array<double, num_equation>& limiting_value) {
#ifdef POST_AI_DATA

	if (std::find(target_cell_indexes_.begin(), target_cell_indexes_.end(), index) == target_cell_indexes_.end())
		return;

	std::string limiting_value_string = "@limiterFunction\n";

	for (size_t i = 0; i < num_equation; ++i)
		limiting_value_string += ms::double_to_str_sp(limiting_value[i]) + "\t";
	limiting_value_string += "\n";

	ai_data_text_set_[index] << std::move(limiting_value_string);

#endif
}




template <size_t space_dimension>
auto PostAI::calculate_face_share_cell_indexes_set(const Grid<space_dimension>& grid) {
	const auto& vnode_index_to_share_cell_indexes = grid.connectivity.vnode_index_to_share_cell_indexes;
	const auto& cell_elements = grid.elements.cell_elements;
	const auto num_cell = cell_elements.size();

	//face share cell indexes set
	std::vector<std::set<size_t>> face_share_cell_indexes_set;
	face_share_cell_indexes_set.reserve(num_cell);

	for (size_t i = 0; i < num_cell; ++i) {
		const auto& element = cell_elements[i];
		const auto& geometry = cell_elements[i].geometry_;

		const auto face_vnode_indexes_set = element.face_vertex_node_indexes_set();
		const auto num_face = face_vnode_indexes_set.size();

		std::set<size_t> face_share_cell_indexes;

		for (const auto& face_vnode_indexes : face_vnode_indexes_set) {
			std::vector<size_t> this_face_share_cell_indexes;

			const auto num_face_vnode = face_vnode_indexes.size();

			const auto& set_0 = vnode_index_to_share_cell_indexes.at(face_vnode_indexes[0]);
			const auto& set_1 = vnode_index_to_share_cell_indexes.at(face_vnode_indexes[1]);
			std::set_intersection(set_0.begin(), set_0.end(), set_1.begin(), set_1.end(), std::back_inserter(this_face_share_cell_indexes));

			if (2 < num_face_vnode) {
				std::vector<size_t> buffer;
				for (size_t i = 2; i < num_face_vnode; ++i) {
					const auto& set_i = vnode_index_to_share_cell_indexes.at(face_vnode_indexes[i]);

					buffer.clear();
					std::set_intersection(this_face_share_cell_indexes.begin(), this_face_share_cell_indexes.end(), set_i.begin(), set_i.end(), std::back_inserter(buffer));
					std::swap(this_face_share_cell_indexes, buffer);
				}
			}

			const auto my_index_pos_iter = std::find(this_face_share_cell_indexes.begin(), this_face_share_cell_indexes.end(), i);
			dynamic_require(my_index_pos_iter != this_face_share_cell_indexes.end(), "my index should be included in this face share cell indexes");

			this_face_share_cell_indexes.erase(my_index_pos_iter);
			dynamic_require(this_face_share_cell_indexes.size() == 1, "face share cell should be unique");

			face_share_cell_indexes.insert(this_face_share_cell_indexes.front());
		}

		face_share_cell_indexes_set.push_back(std::move(face_share_cell_indexes));
	}

	return face_share_cell_indexes_set;
}


template <size_t space_dimension>
auto PostAI::calculate_vertex_nodes_coordinate_string_set(const Grid<space_dimension>& grid) {
	const auto& cell_elements = grid.elements.cell_elements;
	const auto num_cell = cell_elements.size();

	std::vector<std::string> vnodes_coordinate_string_set;
	vnodes_coordinate_string_set.reserve(num_cell);

	for (size_t i = 0; i < num_cell; ++i) {
		const auto& element = cell_elements[i];
		const auto& geometry = element.geometry_;

		const auto vnodes = geometry.vertex_nodes();
		const auto num_vnode = vnodes.size();
		std::string vnodes_coordinate_string = std::to_string(num_vnode) + "\n";

		for (const auto& vnode : vnodes) {
			for (size_t i = 0; i < space_dimension; ++i)
				vnodes_coordinate_string += ms::double_to_str_sp(vnode[i]) + "\t";

			vnodes_coordinate_string += "\n";
		}

		vnodes_coordinate_string_set.push_back(std::move(vnodes_coordinate_string));
	}

	return vnodes_coordinate_string_set;
}


template <size_t num_equation>
auto PostAI::convert_to_solution_strings(const std::vector<EuclideanVector<num_equation>>& solutions) {
	const auto num_solution = solutions.size();
	
	std::vector<std::string> solution_strings;
	solution_strings.reserve(num_solution);



	std::string solution_string;
	for (size_t i = 0; i < num_solution; ++i) {

		const auto& solution = solutions[i];
		for (size_t j = 0; j < num_equation; ++j)
			solution_string += ms::double_to_str_sp(solution[j]) + "\t";

		solution_strings.push_back(std::move(solution_string));
	}

	return solution_strings;
}


