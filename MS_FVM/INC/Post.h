#pragma once
#include "Governing_Equation.h"
#include "Geometry.h"
#include "Text.h"

enum class Post_File_Type {
	Grid, Solution
};


class Post
{
public:
	template <typename Governing_Equation, size_t dim>
	static void intialize(const std::string& name);

	template <size_t dim>
	static void grid(const std::vector<Geometry<dim>>& cell_geometries);

	template <size_t dim>
	static void solution(const std::vector<EuclideanVector<dim>>& solutions, const double time_step, const std::string& comment = "");

private:
	static Text header_text(const Post_File_Type file_type, const double time_step = 0.0);

private:
	static inline std::string name_;
	static inline std::string grid_variables_;
	static inline std::string solution_variables_;
	static inline std::string zone_type_;

	static inline std::vector<size_t> num_post_points_;

	static inline size_t num_element_ = 0;
	static inline size_t num_node_ = 0;	
};


//template definition part
template <typename Governing_Equation, size_t dim>
void Post::intialize(const std::string& name) {
	name_ = name;

	if constexpr (ms::is_Scalar_Eq<Governing_Equation>) {
		if constexpr (dim == 2) {
			grid_variables_ = "Variables = X Y";
			solution_variables_ = "Variables = q";
			zone_type_ = "ZoneType = FETriangle";
		}
	}
	else
		throw std::runtime_error("wrong post initialize");
}

template <size_t dim>
void Post::grid(const std::vector<Geometry<dim>>& cell_geometries) {	
	size_t str_per_line = 1;
	
	const size_t num_cell = cell_geometries.size();
	Post::num_post_points_.resize(num_cell);

	static size_t connectivity_start_index = 1;
		
	Text grid_post_data_text(dim);
	for (size_t i = 0; i < num_cell; ++i) {
		const auto& geometry = cell_geometries[i];

		auto post_nodes = geometry.vertex_nodes();
		Post::num_post_points_[i] = post_nodes.size();
		for (const auto node :post_nodes)
			for (size_t i = 0; i < dim; ++i, ++str_per_line) {
				grid_post_data_text[i] += ms::double_to_string(node[i]) + " ";
				if (str_per_line == 10) {
					grid_post_data_text[i] += "\n";
					str_per_line = 1;
				}

			}

		std::string connectivity_str;
		auto local_connectivities = geometry.local_connectivities();
		for (const auto& local_connectivity : local_connectivities) {
			for (const auto& index : local_connectivity)
				connectivity_str += std::to_string(connectivity_start_index + index) + " ";
			grid_post_data_text << std::move(connectivity_str);
			connectivity_str.clear();	
		}

		connectivity_start_index += Post::num_post_points_[i];
		num_node_ += Post::num_post_points_[i];
		num_element_ += local_connectivities.size();
	}

	auto grid_post_header_text = Post::header_text(Post_File_Type::Grid);

	const auto grid_file_path = "RSC/Post/" + name_ + "/" + name_ + "_grid.plt";
	grid_post_header_text.write(grid_file_path);
	grid_post_data_text.add_write(grid_file_path);
}

template <size_t dim>
void Post::solution(const std::vector<EuclideanVector<dim>>& solutions, const double time_step, const std::string& comment) {
	size_t str_per_line = 1;
	static size_t count = 1;
	
	auto solution_post_header_text = Post::header_text(Post_File_Type::Solution, time_step);

	Text solution_post_data_text(dim);
	for (const auto& solution : solutions){
		for (size_t i = 0; i < Post::num_post_points_[i]; ++i) {
			for (size_t j = 0; j < dim; ++j, ++str_per_line) {
				solution_post_data_text[j] += ms::double_to_string(solution[j]) + " ";
				if (str_per_line == 10) {
					solution_post_data_text[j] += "\n";
					str_per_line = 1;
				}
			}

		}
	}

	std::string solution_file_path;
	if (comment.empty()) 
		solution_file_path = "RSC/Post/" + name_ + "/" + name_ + "_solution_" + std::to_string(count++) + ".plt";
	else
		solution_file_path = "RSC/Post/" + name_ + "/" + name_ + "_" + comment +"_" + std::to_string(count++) + ".plt";

	solution_post_header_text.write(solution_file_path);
	solution_post_data_text.add_write(solution_file_path);
}