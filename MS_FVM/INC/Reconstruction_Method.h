#pragma once
#include <type_traits>
#include <string>

#include "Gradient_Method.h"

class RM {};	// Reconstruction Method


class Constant_Reconstruction : public RM 
{
public:
	static std::string name(void) { return "Constant_Reconstruction"; };
};


template <typename Gradient_Method>
class Linear_Reconstruction : public RM
{
private:
    static constexpr size_t space_dimension_ = Gradient_Method::space_dimension_;

private:
    Gradient_Method gradient_method;

public:
    Linear_Reconstruction(const Grid<space_dimension_>& grid) : gradient_method(grid) {};

    template<size_t num_equation> 
    auto calculate_gradients(const std::vector<EuclideanVector<num_equation>>& solutions) const;

    static std::string name(void) { return "Linear_Reconstruction_" + Gradient_Method::name(); };
};



template <typename Gradient_Method>
class MLP_Base
{
private:
    static constexpr size_t space_dimension_ = Gradient_Method::space_dimension_;

private:
    Gradient_Method gradient_method;

protected:
	std::vector<std::vector<size_t>> vnode_indexes_set_;
	std::vector<Dynamic_Matrix_> center_to_vertex_matrixes_;
	std::unordered_map<size_t, std::set<size_t>> vnode_index_to_share_cell_indexes_;

protected:
	MLP_Base(Grid<space_dimension_>&& grid);

	template<typename Solution>
	std::unordered_map<size_t, std::pair<Solution, Solution>> calculate_vertex_node_index_to_min_max_solution(const std::vector<Solution>& solutions) const;
};


template <typename Gradient_Method>
class MLP_u1 : public MLP_Base<Gradient_Method>
{
private:
    static constexpr size_t space_dimension_ = Gradient_Method::space_dimension_;

public:
    MLP_u1(Grid<space_dimension>&& grid) : MLP_Base<space_dimension>(std::move(grid)) {};

    template<typename Solution>
    void limit_solution_gradients(const std::vector<Solution>& solutions, std::vector<Dynamic_Matrix_>& solution_gradients) const;

    static std::string name(void) { return "MLP_u1_" + Gradient_Method::name(); };


private:
    double limit(const double vertex_solution_delta, const double center_solution, const double min_solution, const double max_solution) const;
};



namespace ms {
	template <typename T>
	inline constexpr bool is_reconsturction_method = std::is_base_of_v<RM, T>;

	template <typename T>
	inline constexpr bool is_MLP_method = std::is_base_of_v<MLP, T>;
}


//template definition part
template <typename Gradient_Method>
template <size_t num_equation>
auto Linear_Reconstruction<Gradient_Method>::calculate_gradients(const std::vector<EuclideanVector<num_equation>>& solutions) const {
    const auto num_cell = solutions.size();

    const auto solution_gradients_temp = gradient_method.calculate_solution_gradients(solutions);

    //dynamic matrix to matrix
    std::vector<Matrix<num_equation, space_dimension_>> solution_gradients;
    solution_gradients.reserve(num_cell);

    for (const auto& solution_gradient : solution_gradients_temp)
        solution_gradients.push_back(solution_gradient);

    return solution_gradients;
}


template <typename Gradient_Method>
MLP_Base<Gradient_Method>::MLP_Base(Grid<space_dimension_>&& grid) : gradient_method(std::move(grid)) {
    SET_TIME_POINT;
    const auto& cell_elements = grid.elements.cell_elements;

    const auto num_cell = cell_elements.size();
    this->vnode_indexes_set_.reserve(num_cell);
    this->center_to_vertex_matrixes_.reserve(num_cell);

    //vnode index to share cell indexes
    this->vnode_index_to_share_cell_indexes_ = std::move(grid.connectivity.vnode_index_to_share_cell_indexes);

    for (size_t i = 0; i < num_cell; ++i) {
        const auto& element = cell_elements[i];
        const auto& geometry = cell_elements[i].geometry_;

        // vnode indexes set
        this->vnode_indexes_set_.push_back(element.vertex_node_indexes());

        //center to vertex matrix
        const auto center_node = geometry.center_node();
        const auto vertex_nodes = geometry.vertex_nodes();
        const auto num_vertex = vertex_nodes.size();

        Dynamic_Matrix_ center_to_vertex_matrix(space_dimension_, num_vertex);
        for (size_t i = 0; i < num_vertex; ++i) {
            const auto center_to_vertex = vertex_nodes[i] - center_node;
            for (size_t j = 0; j < space_dimension_; ++j)
                center_to_vertex_matrix.at(j, i) = center_to_vertex[j];
        }
        this->center_to_vertex_matrixes_.push_back(std::move(center_to_vertex_matrix));
    }

    Log::content_ << std::left << std::setw(50) << "@ Construct Cells FVM MLP Base" << " ----------- " << GET_TIME_DURATION << "s\n\n";
    Log::print();
}

template <typename Gradient_Method>
template <typename Solution>
std::unordered_map<size_t, std::pair<Solution, Solution>> MLP_Base<Gradient_Method>::calculate_vertex_node_index_to_min_max_solution(const std::vector<Solution>& solutions) const {
    constexpr size_t num_equation = Solution::dimension();
    const size_t num_vnode = this->vnode_index_to_share_cell_indexes_.size();

    std::unordered_map<size_t, std::pair<Solution, Solution>> vnode_index_to_min_max_solution;
    vnode_index_to_min_max_solution.reserve(num_vnode);

    for (const auto& [vnode_index, share_cell_indexes] : this->vnode_index_to_share_cell_indexes_) {
        const size_t num_share_cell = share_cell_indexes.size();
        std::array<std::vector<double>, num_equation> equation_wise_solutions;

        for (size_t i = 0; i < num_equation; ++i)
            equation_wise_solutions[i].reserve(num_share_cell);

        for (const auto cell_index : share_cell_indexes) {
            for (size_t i = 0; i < num_equation; ++i)
                equation_wise_solutions[i].push_back(solutions[cell_index][i]);
        }

        std::array<double, num_equation> min_solution;
        std::array<double, num_equation> max_solution;
        for (size_t i = 0; i < num_equation; ++i) {
            min_solution[i] = *std::min_element(equation_wise_solutions[i].begin(), equation_wise_solutions[i].end());
            max_solution[i] = *std::max_element(equation_wise_solutions[i].begin(), equation_wise_solutions[i].end());
        }
        Solution min_sol = min_solution;
        Solution max_sol = max_solution;

        vnode_index_to_min_max_solution.emplace(vnode_index, std::make_pair(min_sol, max_sol));
    }

    return vnode_index_to_min_max_solution;
}

template <typename Gradient_Method>
template <typename Solution>
void MLP_u1<Gradient_Method>::limit_solution_gradients(const std::vector<Solution>& solutions, std::vector<Dynamic_Matrix_>& solution_gradients) const {
    constexpr size_t num_equation = Solution::dimension();
    const auto num_cell = solutions.size();

    const auto vnode_index_to_min_max_solution = this->calculate_vertex_node_index_to_min_max_solution(solutions);

    for (size_t i = 0; i < num_cell; ++i) {
        auto& gradient = solution_gradients[i];
        const auto vertex_solution_delta_matrix = gradient * this->center_to_vertex_matrixes_[i];

        std::array<double, num_equation> limiting_values;
        limiting_values.fill(1);

        const auto& vnode_indexes = this->vnode_indexes_set_[i];
        const auto num_vertex = vnode_indexes.size();

        for (size_t j = 0; j < num_vertex; ++j) {
            const auto vnode_index = vnode_indexes[j];
            const auto& [min_solution, max_solution] = vnode_index_to_min_max_solution.at(vnode_index);

            for (size_t e = 0; e < num_equation; ++e) {
                const auto limiting_value = this->limit(vertex_solution_delta_matrix.at(e, j), solutions[i].at(e), min_solution.at(e), max_solution.at(e));
                limiting_values[e] = min(limiting_values[e], limiting_value);
            }
        }

        for (size_t i = 0; i < num_equation; ++i)
            for (size_t j = 0; j < space_dimension; ++j)
                gradient.at(i, j) *= limiting_values.at(i);
    }
}

template <typename Gradient_Method>
double MLP_u1<Gradient_Method>::limit(const double vertex_solution_delta, const double center_solution, const double min_solution, const double max_solution) const {
    if (vertex_solution_delta < 0)
        return min((min_solution - center_solution) / vertex_solution_delta, 1);
    else
        return min((max_solution - center_solution) / vertex_solution_delta, 1);
}