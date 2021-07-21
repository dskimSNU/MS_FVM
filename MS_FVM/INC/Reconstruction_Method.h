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


template<size_t num_equation, size_t space_dimension>
struct Linear_Reconstructed_Solution
{
private:
    using Solution_             = EuclideanVector<num_equation>;
    using Solution_Gradient_    = Matrix<num_equation, space_dimension>;

public:
    const std::vector<Solution_>&   solutions;
    std::vector<Solution_Gradient_> solution_gradients;
};


template <typename Gradient_Method>
class Linear_Reconstruction : public RM
{
private:
    static constexpr size_t num_equation_       = Gradient_Method::num_equation_;
    static constexpr size_t space_dimension_    = Gradient_Method::space_dimension_;

private:
    Gradient_Method gradient_method;

public:
    Linear_Reconstruction(const Grid<space_dimension_>& grid) : gradient_method(grid) {};


    auto reconstruct_solutions(const std::vector<EuclideanVector<num_equation_>>& solutions) const;

    static std::string name(void) { return "Linear_Reconstruction_" + Gradient_Method::name(); };
};



template <typename Gradient_Method>
class MLP_Base : public RM
{
private:
    static constexpr size_t num_equation_       = Gradient_Method::num_equation_;
    static constexpr size_t space_dimension_    = Gradient_Method::space_dimension_;

    using Solution_ = EuclideanVector<num_equation_>;

protected:
    Gradient_Method gradient_method;

	std::vector<std::vector<size_t>> vnode_indexes_set_;
	std::vector<Dynamic_Matrix_> center_to_vertex_matrixes_;
	std::unordered_map<size_t, std::set<size_t>> vnode_index_to_share_cell_indexes_;

public:
    auto reconstruct_solutions(const std::vector<Solution_>& solutions) const;

protected:
    MLP_Base(Grid<space_dimension_>&& grid);

	auto calculate_vertex_node_index_to_min_max_solution(const std::vector<Solution_>& solutions) const;

    virtual double limit(const double vertex_solution_delta, const double center_solution, const double min_solution, const double max_solution) const abstract;
};


template <typename Gradient_Method>
class MLP_u1 : public MLP_Base<Gradient_Method>
{
    static constexpr size_t space_dimension_ = Gradient_Method::space_dimension_;

public:
    MLP_u1(Grid<space_dimension_>&& grid) : MLP_Base<Gradient_Method>(std::move(grid)) {};

    static std::string name(void) { return "MLP_u1_" + Gradient_Method::name(); };

private:
    double limit(const double vertex_solution_delta, const double center_solution, const double min_solution, const double max_solution) const override;
};



namespace ms {
	template <typename T>
	inline constexpr bool is_reconsturction_method = std::is_base_of_v<RM, T>;
}


//template definition part
template <typename Gradient_Method>
auto Linear_Reconstruction<Gradient_Method>::reconstruct_solutions(const std::vector<EuclideanVector<num_equation_>>& solutions) const {
    const auto num_cell = solutions.size();
    const auto solution_gradients_temp = gradient_method.calculate_solution_gradients(solutions);

    //dynamic matrix to matrix
    std::vector<Matrix<num_equation_, space_dimension_>> solution_gradients;
    solution_gradients.reserve(num_cell);

    for (const auto& solution_gradient : solution_gradients_temp)
        solution_gradients.push_back(solution_gradient);

    return Linear_Reconstructed_Solution<num_equation_, space_dimension_>{ solutions, solution_gradients };
}


template <typename Gradient_Method>
auto MLP_Base<Gradient_Method>::reconstruct_solutions(const std::vector<Solution_>& solutions) const {
    auto solution_gradients = this->gradient_method.calculate_solution_gradients(solutions);
    const auto vnode_index_to_min_max_solution = this->calculate_vertex_node_index_to_min_max_solution(solutions);

    const auto num_cell = solutions.size();
    for (size_t i = 0; i < num_cell; ++i) {
        auto& gradient = solution_gradients[i];
        const auto vertex_solution_delta_matrix = gradient * this->center_to_vertex_matrixes_[i];

        std::array<double, num_equation_> limiting_values;
        limiting_values.fill(1);

        const auto& vnode_indexes = this->vnode_indexes_set_[i];
        const auto num_vertex = vnode_indexes.size();

        for (size_t j = 0; j < num_vertex; ++j) {
            const auto vnode_index = vnode_indexes[j];
            const auto& [min_solution, max_solution] = vnode_index_to_min_max_solution.at(vnode_index);

            for (size_t e = 0; e < num_equation_; ++e) {
                const auto limiting_value = this->limit(vertex_solution_delta_matrix.at(e, j), solutions[i].at(e), min_solution.at(e), max_solution.at(e));
                limiting_values[e] = min(limiting_values[e], limiting_value);
            }
        }

        for (size_t i = 0; i < num_equation_; ++i)
            for (size_t j = 0; j < space_dimension_; ++j)
                gradient.at(i, j) *= limiting_values.at(i);
    }

    //dynamic matrix to matrix
    std::vector<Matrix<num_equation_, space_dimension_>> limited_solution_gradient;
    limited_solution_gradient.reserve(num_cell);

    for (const auto& solution_gradient : solution_gradients)
        limited_solution_gradient.push_back(solution_gradient);

    return Linear_Reconstructed_Solution<num_equation_, space_dimension_>{ solutions, limited_solution_gradient };
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
auto MLP_Base<Gradient_Method>::calculate_vertex_node_index_to_min_max_solution(const std::vector<Solution_>& solutions) const {
    const size_t num_vnode = this->vnode_index_to_share_cell_indexes_.size();

    std::unordered_map<size_t, std::pair<Solution_, Solution_>> vnode_index_to_min_max_solution;
    vnode_index_to_min_max_solution.reserve(num_vnode);

    for (const auto& [vnode_index, share_cell_indexes] : this->vnode_index_to_share_cell_indexes_) {
        const size_t num_share_cell = share_cell_indexes.size();
        std::array<std::vector<double>, num_equation_> equation_wise_solutions;

        for (size_t i = 0; i < num_equation_; ++i)
            equation_wise_solutions[i].reserve(num_share_cell);

        for (const auto cell_index : share_cell_indexes) {
            for (size_t i = 0; i < num_equation_; ++i)
                equation_wise_solutions[i].push_back(solutions[cell_index][i]);
        }

        std::array<double, num_equation_> min_solution;
        std::array<double, num_equation_> max_solution;
        for (size_t i = 0; i < num_equation_; ++i) {
            min_solution[i] = *std::min_element(equation_wise_solutions[i].begin(), equation_wise_solutions[i].end());
            max_solution[i] = *std::max_element(equation_wise_solutions[i].begin(), equation_wise_solutions[i].end());
        }
        Solution_ min_sol = min_solution;
        Solution_ max_sol = max_solution;

        vnode_index_to_min_max_solution.emplace(vnode_index, std::make_pair(min_sol, max_sol));
    }

    return vnode_index_to_min_max_solution;
}


template <typename Gradient_Method>
double MLP_u1<Gradient_Method>::limit(const double vertex_solution_delta, const double center_solution, const double min_solution, const double max_solution) const {
    if (vertex_solution_delta < 0)
        return min((min_solution - center_solution) / vertex_solution_delta, 1);
    else
        return min((max_solution - center_solution) / vertex_solution_delta, 1);
}