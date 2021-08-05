#pragma once
#include <type_traits>
#include <string>

#include "Gradient_Method.h"
#include "PostAI.h"

class RM {};	// Reconstruction Method


class Constant_Reconstruction : public RM 
{
public:
    template <size_t space_dimension>
    Constant_Reconstruction(const Grid<space_dimension>& grid) {}; //because semi discrete equation constructor

    static std::string name(void) { return "Constant_Reconstruction"; };
};


template<size_t num_equation, size_t space_dimension>
struct Linear_Reconstructed_Solution
{
private:
    using Solution_             = Euclidean_Vector<num_equation>;
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


    auto reconstruct_solutions(const std::vector<Euclidean_Vector<num_equation_>>& solutions) const;

    static std::string name(void) { return "Linear_Reconstruction_" + Gradient_Method::name(); };
};



template <typename Gradient_Method>
class MLP_Base : public RM
{
private:
    static constexpr size_t num_equation_       = Gradient_Method::num_equation_;
    static constexpr size_t space_dimension_    = Gradient_Method::space_dimension_;

    using Solution_ = Euclidean_Vector<num_equation_>;

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


template <typename Gradient_Method>
class ANN_limiter : public RM
{
private:
    static constexpr size_t num_equation_ = Gradient_Method::num_equation_;
    static constexpr size_t space_dimension_ = Gradient_Method::space_dimension_;

    using Solution_ = Euclidean_Vector<num_equation_>;

//private: for test
public:
    Gradient_Method gradient_method_;

    std::vector<std::set<size_t>> vertex_share_cell_index_sets;
    std::vector<std::set<size_t>> face_share_cell_index_sets;

public:
    ANN_limiter(const Grid<space_dimension_>& grid);

    auto reconstruct_solutions(const std::vector<Solution_>& solutions);
    static std::string name(void) { return "ANN_Reconstruction_" + Gradient_Method::name(); };

private:
    auto calculate_set_of_face_share_cell_index_set(const Grid<space_dimension_>& grid) const;
    std::vector<size_t> ordering_function(const std::vector<Solution_>& solutions, const size_t target_cell_index) const;
    bool is_constant_region(const std::vector<Solution_>& solutions, const size_t target_cell_index) const;
    double limit(const double* ptr) const;
};



namespace ms {
	template <typename T>
	inline constexpr bool is_reconsturction_method = std::is_base_of_v<RM, T>;

    template <typename T>
    inline constexpr bool is_constant_reconstruction = std::is_same_v<Constant_Reconstruction, T>;
}


//template definition part
template <typename Gradient_Method>
auto Linear_Reconstruction<Gradient_Method>::reconstruct_solutions(const std::vector<Euclidean_Vector<num_equation_>>& solutions) const {
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

    PostAI::record_solution_datas(solutions, solution_gradients);

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

        PostAI::record_limiting_value(i, limiting_values);

        for (size_t i = 0; i < num_equation_; ++i)
            for (size_t j = 0; j < space_dimension_; ++j)
                gradient.at(i, j) *= limiting_values.at(i);
    }

    PostAI::post();

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

template <typename Gradient_Method>
ANN_limiter<Gradient_Method>::ANN_limiter(const Grid<space_dimension_>& grid) : gradient_method_(grid) {
    this->face_share_cell_index_sets = this->calculate_set_of_face_share_cell_index_set(grid);

    const auto& vnode_index_to_share_cell_indexes = grid.connectivity.vnode_index_to_share_cell_indexes;
    const auto& cell_elements = grid.elements.cell_elements;

    const auto num_cell = cell_elements.size();
    this->vertex_share_cell_index_sets.reserve(num_cell);    

    for (size_t i = 0; i < num_cell; ++i) {
        const auto& cell_element = cell_elements[i];
        const auto& cell_geometry = cell_element.geometry_;

        
        std::set<size_t> vertex_share_cell_index_set;

        const auto vnode_indexes = cell_element.vertex_node_indexes();
        for (const auto& vnode_index : vnode_indexes) {
            const auto& vnode_share_cell_index_set = vnode_index_to_share_cell_indexes.at(vnode_index);
            vertex_share_cell_index_set.insert(vnode_share_cell_index_set.begin(), vnode_share_cell_index_set.end());
        }

        this->face_share_cell_index_sets.push_back(std::move(vertex_share_cell_index_set));
    }
}

template <typename Gradient_Method>
auto ANN_limiter<Gradient_Method>::reconstruct_solutions(const std::vector<Solution_>& solutions) {
    static const auto num_solution = solutions.size();
    
    auto solution_gradients = this->gradient_method_.calculate_solution_gradients(solutions);

    for (size_t i = 0; i < num_solution; ++i) {
        if (this->is_constant_region(solutions, i))
            continue;
                
        const auto ordered_indexes = this->ordering_function(solutions, i);

        const auto num_vertex_share_cell = this->vertex_share_cell_index_sets.at(i).size();
        std::vector<double> input(num_vertex_share_cell * 3);

        for (size_t j = 0; j < num_equation_; ++j) {
            for (size_t k = 0; k < num_vertex_share_cell; ++k) {
                const auto index = ordered_indexes[k];

                const auto& solution = solutions[index];
                const auto& solution_gradient = solution_gradients[index];

                const auto solution_start_index = 0;
                const auto solution_gradient_x_start_index = num_vertex_share_cell;
                const auto solution_gradient_y_start_index = 2 * num_vertex_share_cell;

                input[solution_start_index + k] = solution.at(j);
                input[solution_gradient_x_start_index + k] = solution_gradient.at(j, 0);
                input[solution_gradient_y_start_index + k] = solution_gradient.at(j, 1);
            }
        }

        auto limiter_value = this->limit(input.data());

        auto& solution_gradient = solution_gradients[i];

        for (size_t i = 0; i < num_equation_; ++i)
            for (size_t j = 0; j < space_dimension_; ++j)
                solution_gradient.at(i, j) *= limiter_value;
    }

    //dynamic matrix to matrix
    std::vector<Matrix<num_equation_, space_dimension_>> limited_solution_gradient;
    limited_solution_gradient.reserve(num_solution);

    for (const auto& solution_gradient : solution_gradients)
        limited_solution_gradient.push_back(solution_gradient);

    return Linear_Reconstructed_Solution<num_equation_, space_dimension_>{ solutions, limited_solution_gradient };
}

template <typename Gradient_Method>
auto ANN_limiter<Gradient_Method>::calculate_set_of_face_share_cell_index_set(const Grid<space_dimension_>& grid) const {
    const auto& vnode_index_to_share_cell_indexes = grid.connectivity.vnode_index_to_share_cell_indexes;
    const auto& cell_elements = grid.elements.cell_elements;
    const auto num_cell = cell_elements.size();

    //face share cell indexes set
    std::vector<std::set<size_t>> set_of_face_share_cell_index_set;
    set_of_face_share_cell_index_set.reserve(num_cell);

    for (size_t i = 0; i < num_cell; ++i) {
        const auto& element = cell_elements[i];
        const auto& geometry = cell_elements[i].geometry_;

        const auto set_of_face_vnode_indexes = element.face_vertex_node_indexes_set();
        const auto num_face = set_of_face_vnode_indexes.size();

        std::set<size_t> face_share_cell_index_set;

        for (const auto& face_vnode_indexes : set_of_face_vnode_indexes) {
            std::vector<size_t> this_face_share_cell_indexes;

            const auto& set_0 = vnode_index_to_share_cell_indexes.at(face_vnode_indexes[0]);
            const auto& set_1 = vnode_index_to_share_cell_indexes.at(face_vnode_indexes[1]);
            std::set_intersection(set_0.begin(), set_0.end(), set_1.begin(), set_1.end(), std::back_inserter(this_face_share_cell_indexes));


            const auto num_face_vnode = face_vnode_indexes.size();
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

            face_share_cell_index_set.insert(this_face_share_cell_indexes.front());
        }

        set_of_face_share_cell_index_set.push_back(std::move(face_share_cell_index_set));
    }

    return set_of_face_share_cell_index_set;
}


template <typename Gradient_Method>
std::vector<size_t> ANN_limiter<Gradient_Method>::ordering_function(const std::vector<Solution_>& solutions, const size_t target_cell_index) const {
    const auto& vnode_share_cell_index_set = this->vertex_share_cell_index_sets.at(target_cell_index);

    const auto num_cell = vnode_share_cell_index_set.size();
    std::vector<size_t> ordered_indexes;
    ordered_indexes.reserve(num_cell);

    ordered_indexes.push_back(target_cell_index);

    while (ordered_indexes.size() == num_cell) {
        const auto& face_share_cell_index_set = this->face_share_cell_index_sets.at(ordered_indexes.back());

        std::vector<size_t> face_share_cell_indexes_in_chunk;
        std::set_intersection(vnode_share_cell_index_set.begin(), vnode_share_cell_index_set.end(), face_share_cell_index_set.begin(), face_share_cell_index_set.end(), std::back_inserter(face_share_cell_indexes_in_chunk));

        std::vector<size_t> candidate_cell_indexes;
        std::set<size_t> ordered_index_set(ordered_indexes.begin(), ordered_indexes.end());
        std::set_difference(face_share_cell_indexes_in_chunk.begin(), face_share_cell_indexes_in_chunk.end(), ordered_index_set.begin(), ordered_index_set.end(), std::back_inserter(candidate_cell_indexes));

        if (1 < candidate_cell_indexes.size()) {
            std::vector<double> temporary_solutions;
            temporary_solutions.reserve(candidate_cell_indexes.size());

            for (const auto candidate_cell_index : candidate_cell_indexes)
                temporary_solutions.push_back(solutions[candidate_cell_index].at(0));

            const auto max_solution_iter = std::max_element(temporary_solutions.begin(), temporary_solutions.end());
            const auto pos = max_solution_iter - temporary_solutions.begin();

            const auto index = *std::next(candidate_cell_indexes.begin(), pos);
            ordered_indexes.push_back(index);
        }
        else
            ordered_indexes.push_back(candidate_cell_indexes.front());
    }

    return ordered_indexes;
}

template <typename Gradient_Method>
bool ANN_limiter<Gradient_Method>::is_constant_region(const std::vector<Solution_>& solutions, const size_t target_cell_index) const {
    const auto& vertex_share_cell_index_set = this->vertex_share_cell_index_sets.at(target_cell_index);
    const auto num_vertex_share_cell = vertex_share_cell_index_set.size();

    std::vector<double> vertex_share_cell_solutions;
    vertex_share_cell_solutions.reserve(num_vertex_share_cell);

    for (const auto vertex_share_cell_index : vertex_share_cell_index_set)
        vertex_share_cell_solutions.push_back(solutions[vertex_share_cell_index][0]);

    const auto min_solution = *std::min_element(vertex_share_cell_solutions.begin(), vertex_share_cell_solutions.end());
    const auto max_solution = *std::max_element(vertex_share_cell_solutions.begin(), vertex_share_cell_solutions.end());
    const auto solution_diff = max_solution - min_solution;

    if (solution_diff < 0.01)
        return true;
    else
        return false;
}