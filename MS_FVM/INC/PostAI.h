#pragma once
#include "Grid_Builder.h"
#include "Governing_Equation.h"
#include "Inital_Condition.h"

class PostAI2D
{
private:
	inline static double x_advection_speed_;
	inline static double y_advection_speed_;
	inline static std::vector<std::set<size_t>> chunk_cell_indexes_set_;
	inline static std::vector<std::set<std::set<size_t>>> face_connectivities_set_;
	inline static std::vector<std::vector<EuclideanVector<2>>> vnodes_set_;
	inline static std::vector<EuclideanVector<2>> cnodes_;
	inline static const double* time_ptr_;
	inline static std::unordered_set<size_t> target_cell_indexes_;

public:
	static void calculate_target_cell_indexes(void) {
		const auto t = *time_ptr_;

		//Assume that domian [0,1] x [0,1]
		const auto exact_x_start	= 0.25 + x_advection_speed_ * t - static_cast<int>(0.25 + x_advection_speed_ * t);
		const auto exact_x_end		= 0.75 + x_advection_speed_ * t - static_cast<int>(0.75 + x_advection_speed_ * t);
		const auto exact_y_start	= 0.25 + y_advection_speed_ * t - static_cast<int>(0.25 + y_advection_speed_ * t);
		const auto exact_y_end		= 0.75 + y_advection_speed_ * t - static_cast<int>(0.75 + y_advection_speed_ * t);



	}

	static void syncrosize_time(const double& current_time) {
		time_ptr_ = &current_time;
	}

	template <typename Governing_Equation, typename Initial_Condition>
	static void intialize(const Grid<2>& grid) {
		static_require(std::is_same_v<Governing_Equation, Linear_Advection_2D>, "Governing Equation should be Linear advection 2D");
		static_require(std::is_same_v<Initial_Condition, Square_Wave_2D>,		"Initial Colndition should be Square wave 2D");


		const auto [x_advection_speed, y_advection_speed] = Linear_Advection_2D::advection_speed();
		x_advection_speed_ = x_advection_speed;
		y_advection_speed_ = y_advection_speed;


		const auto& vnode_index_to_share_cell_indexes = grid.connectivity.vnode_index_to_share_cell_indexes;
		const auto& cell_elements = grid.elements.cell_elements;
		const auto num_cell = cell_elements.size();


		//face share cell indexes set
		std::vector<std::vector<size_t>> face_share_cell_indexes_set;
		face_share_cell_indexes_set.reserve(num_cell);

		for (size_t i = 0; i < num_cell; ++i) {
			const auto& element = cell_elements[i];
			const auto& geometry = cell_elements[i].geometry_;

			const auto face_vnode_indexes_set = element.face_vertex_node_indexes_set();
			const auto num_face = face_vnode_indexes_set.size();

			std::vector<size_t> face_share_cell_indexes;
			face_share_cell_indexes.reserve(num_face);

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

				face_share_cell_indexes.push_back(this_face_share_cell_indexes.front());
			}

			face_share_cell_indexes_set.push_back(std::move(face_share_cell_indexes));
		}




		chunk_cell_indexes_set_.reserve(num_cell);
		face_connectivities_set_.reserve(num_cell);
		vnodes_set_.reserve(num_cell);
		cnodes_.reserve(num_cell);

		for (size_t i = 0; i < num_cell; ++i) {
			const auto& cell_element = cell_elements[i];
			const auto& cell_geometry = cell_element.geometry_;

			//chunk_cell_indexes_set_
			std::set<size_t> chunk_cell_indexes;
			
			const auto vnode_indexes = cell_element.vertex_node_indexes();
			for (const auto& vnode_index : vnode_indexes) {
				const auto& vnode_share_cell_indexes = vnode_index_to_share_cell_indexes.at(vnode_index);
				chunk_cell_indexes.insert(vnode_share_cell_indexes.begin(), vnode_share_cell_indexes.end());
			}
			chunk_cell_indexes.erase(i);
			
			//face_connectivities_set_
			std::set<std::set<size_t>> face_connetivities;

			for (const auto chunk_cell_index : chunk_cell_indexes) {
				const auto& face_share_cell_indexes = face_share_cell_indexes_set.at(chunk_cell_index); //sorted vector
				
				std::vector<size_t> face_share_cell_indexes_in_chunk;
				std::set_intersection(chunk_cell_indexes.begin(), chunk_cell_indexes.end(), face_share_cell_indexes.begin(), face_share_cell_indexes.end(), std::back_inserter(face_share_cell_indexes_in_chunk));

				for (const auto face_share_cell_index_in_chunk : face_share_cell_indexes_in_chunk)
					face_connetivities.insert({ chunk_cell_index, face_share_cell_index_in_chunk });
			}

			//vnodes_set
			const auto vnodes = cell_geometry.vertex_nodes();
			//const auto num_node = vnodes.size();

			//	//container change
			//std::vector<std::vector<double>> vnodes_vec(num_node);
			//for (size_t i = 0; i < num_node; ++i) {
			//	const auto& vnode = vnodes[i];
			//	auto& vnode_vec = vnodes_vec[i];
			//	vnode_vec.resize(space_dimension);

			//	for (size_t j = 0; j < space_dimension; ++j)
			//		vnode_vec[j] = vnode[j];
			//}

			//cell_centers
			const auto cnode = cell_geometry.center_node();

			//	//container change
			//std::vector<double> cnode_vec(space_dimension);
			//for (size_t i = 0; i < space_dimension; ++i)
			//	cnode_vec[i] = cnode[i];


			chunk_cell_indexes_set_.push_back(std::move(chunk_cell_indexes));
			face_connectivities_set_.push_back(std::move(face_connetivities));
			vnodes_set_.push_back(std::move(vnodes));
			cnodes_.push_back(std::move(cnode));

			//vnodes_set_.push_back(std::move(vnodes_vec));
			//cnodes_.push_back(std::move(cnode_vec));
		}
	}	
};