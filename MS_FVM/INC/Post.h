#pragma once
#include "Governing_Equation.h"
#include "Element.h"
#include "Text.h"

enum class Post_File_Type {
	Grid, Solution
};

template <typename Governing_Equation>
class Post
{
private:
	static_require(ms::is_governing_equation<Governing_Equation>, "it should be governing equation");

	static constexpr size_t space_dimension_	= Governing_Equation::space_dimension();
	static constexpr size_t num_equation_		= Governing_Equation::num_equation();

private:
	static inline std::string path_;
	static inline std::string grid_variable_str_;
	static inline std::string solution_variable_str_;
	static inline std::string zone_type_str_;
	static inline std::vector<size_t> num_post_points_;

	static inline size_t num_element_ = 0;
	static inline size_t num_node_ = 0;

	static inline const double* time_ptr_ = nullptr;
public:
	static void set_path(const std::string& path) { Post::path_ = path; };	
	static void intialize(void);	
	static void grid(const std::vector<Element<space_dimension_>>& cell_elements);	
	static void syncronize_time(const double& current_time) { Post::time_ptr_ = &current_time; };
	static void solution(const std::vector<EuclideanVector<num_equation_>>& solutions, const std::string& comment = "");

private:
	static Text header_text(const Post_File_Type file_type);
};


//template definition part
template <typename Governing_Equation>
void Post<Governing_Equation>::intialize(void) {
	if constexpr (ms::is_SCL_2D<Governing_Equation>) {
			grid_variable_str_ = "Variables = X Y";
			solution_variable_str_ = "Variables = q";
			zone_type_str_ = "ZoneType = FETriangle";
	}
	else if constexpr (std::is_same_v<Governing_Equation, Euler_2D>) {
		grid_variable_str_ = "Variables = X Y";
		solution_variable_str_ = "Variables = rho rhou rhov rhoE u v p a";
		zone_type_str_ = "ZoneType = FETriangle";
	}
	else 
		throw std::runtime_error("wrong post initialize");
}

template <typename Governing_Equation>
void Post<Governing_Equation>::grid(const std::vector<Element<space_dimension_>>& cell_elements) {
	const size_t num_cell = cell_elements.size();
	Post::num_post_points_.resize(num_cell);

	size_t str_per_line = 1;
	size_t connectivity_start_index = 1;

	Text grid_post_data_text(space_dimension_);
	for (size_t i = 0; i < num_cell; ++i) {
		const auto& geometry = cell_elements[i].geometry_;

		auto post_nodes = geometry.vertex_nodes();
		Post::num_post_points_[i] = post_nodes.size();
		for (const auto node : post_nodes) {
			for (size_t i = 0; i < space_dimension_; ++i, ++str_per_line) {
				grid_post_data_text[i] += ms::double_to_string(node[i]) + " ";
				if (str_per_line == 10) {
					grid_post_data_text[i] += "\n";
					str_per_line = 1;
				}
			}
		}

		std::string connectivity_str;
		auto local_connectivities = geometry.reference_geometry_.local_connectivities();
		for (const auto& local_connectivity : local_connectivities) {
			for (const auto& index : local_connectivity)
				connectivity_str += std::to_string(connectivity_start_index + index) + " ";
			grid_post_data_text << std::move(connectivity_str);
		}

		connectivity_start_index += Post::num_post_points_[i];
		num_node_ += Post::num_post_points_[i];
		num_element_ += local_connectivities.size();
	}

	auto grid_post_header_text = Post::header_text(Post_File_Type::Grid);

	const auto grid_file_path = path_ + "grid.plt";
	grid_post_header_text.write(grid_file_path);
	grid_post_data_text.add_write(grid_file_path);
}

template <typename Governing_Equation>
void Post<Governing_Equation>::solution(const std::vector<EuclideanVector<num_equation_>>& solutions, const std::string& comment) {
	static size_t count = 1;

	std::string solution_file_path;
	if (comment.empty())
		solution_file_path = path_ + "solution_" + std::to_string(count++) + ".plt";
	else
		solution_file_path = path_ + "solution_" + std::to_string(count++) + "_" + comment + ".plt";

	//solution post header text
	auto solution_post_header_text = Post::header_text(Post_File_Type::Solution);
	solution_post_header_text.write(solution_file_path);

	
	//solution post data text
	size_t str_per_line = 1;	

	Text solution_post_data_text;

	if constexpr (ms::is_SCL_2D<Governing_Equation>) {
		solution_post_data_text.resize(num_equation_);
		
		const auto num_solution = solutions.size();
		for (size_t i = 0; i < num_solution; ++i) {
			const auto& solution = solutions[i];
			for (size_t j = 0; j < Post::num_post_points_[i]; ++j, ++str_per_line) {
				solution_post_data_text[0] += ms::double_to_string(solution[0]) + " ";
				if (str_per_line == 10) {
					solution_post_data_text[0] += "\n";
					str_per_line = 1;
				}
			}
		}
	}
	else if constexpr (std::is_same_v<Governing_Equation, Euler_2D>) {
		solution_post_data_text.resize(2 * num_equation_);

		const auto num_solution = solutions.size();
		for (size_t i = 0; i < num_solution; ++i) {
			const auto& cvariable = solutions[i];
			const auto pvariable = Euler_2D::conservative_to_primitive(cvariable);

			for (size_t j = 0; j < Post::num_post_points_[i]; ++j) {
				for (size_t k = 0; k < num_equation_; ++k, ++str_per_line) {
					solution_post_data_text[k] += ms::double_to_string(cvariable[k]) + " ";
					if (str_per_line == 10) {
						solution_post_data_text[k] += "\n";
						str_per_line = 1;
					}
				}

				for (size_t k = 0; k < num_equation_; ++k, ++str_per_line) {
					solution_post_data_text[k + 4] += ms::double_to_string(pvariable[k]) + " ";
					if (str_per_line == 10) {
						solution_post_data_text[k + 4] += "\n";
						str_per_line = 1;
					}
				}
			}
		}
	}

	solution_post_data_text.add_write(solution_file_path);
}

template <typename Governing_Equation>
Text Post<Governing_Equation>::header_text(const Post_File_Type file_type) {
	static size_t strand_id = 0;

	Text header;
	header.reserve(10);
	if (file_type == Post_File_Type::Grid) {
		header << "Title = Grid";
		header << "FileType = Grid";
		header << grid_variable_str_;
		header << "Zone T = Grid";
	}
	else {
		header << "Title = Solution_at_" + ms::double_to_string(*time_ptr_);
		header << "FileType = Solution";
		header << solution_variable_str_;
		strand_id++;
		header << "Zone T = Solution_at_" + ms::double_to_string(*time_ptr_);
	}

	header << zone_type_str_;
	header << "Nodes = " + std::to_string(num_node_);
	header << "Elements = " + std::to_string(num_element_);
	header << "DataPacking = Block";
	header << "StrandID = " + std::to_string(strand_id);

	if (file_type == Post_File_Type::Grid)
		header << "SolutionTime = 0.0 \n\n";
	else
		header << "SolutionTime = " + ms::double_to_string(*time_ptr_) + "\n\n";
	
	return header;
}