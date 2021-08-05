#pragma once
#include "Grid_Builder.h"

template <size_t num_equation, size_t space_dimension>
class Least_Square_Base
{
public:    
    static constexpr size_t space_dimension_    = space_dimension;
    static constexpr size_t num_equation_       = num_equation;
    
    using Solution_ = Euclidean_Vector<num_equation>;

protected:
    size_t num_cell_;
    std::vector<std::vector<size_t>> near_cell_indexes_set_;
    std::vector<Dynamic_Matrix_> least_square_matrixes_;

public:
    std::vector<Dynamic_Matrix_> calculate_solution_gradients(const std::vector<Solution_>& solutions) const;

protected:
    std::vector<Dynamic_Matrix_> calculate_solution_delta_matrixes(const std::vector<Solution_>& solutions) const;
};


template <size_t num_equation, size_t space_dimension>
class Vertex_Least_Square : public Least_Square_Base<num_equation, space_dimension>
{
public:
    Vertex_Least_Square(const Grid<space_dimension>& grid);

    static std::string name(void) { return "Vertex_Least_Square"; };
};


template <size_t num_equation, size_t space_dimension>
class Face_Least_Square : public Least_Square_Base<num_equation, space_dimension>
{
public:
    Face_Least_Square(const Grid<space_dimension>& grid);

    static std::string name(void) { return "Face_Least_Square"; };
};



//template definition part
template <size_t num_equation, size_t space_dimension>
std::vector<Dynamic_Matrix_> Least_Square_Base<num_equation, space_dimension>::calculate_solution_gradients(const std::vector<Solution_>& solutions) const {
    std::vector<Dynamic_Matrix_> solution_gradients;
    solution_gradients.reserve(this->num_cell_);

    const auto solution_delta_matrixes = this->calculate_solution_delta_matrixes(solutions);

    for (size_t i = 0; i < this->num_cell_; ++i)
        solution_gradients.push_back(solution_delta_matrixes[i] * this->least_square_matrixes_[i]);

    return solution_gradients;
}

template <size_t num_equation, size_t space_dimension>
std::vector<Dynamic_Matrix_> Least_Square_Base<num_equation, space_dimension>::calculate_solution_delta_matrixes(const std::vector<Solution_>& solutions) const {
    std::vector<Dynamic_Matrix_> solution_delta_matrixes;
    solution_delta_matrixes.reserve(this->num_cell_);

    for (size_t i = 0; i < this->num_cell_; ++i) {
        const auto& near_cell_indexes = this->near_cell_indexes_set_.at(i);
        const auto num_near_cell = near_cell_indexes.size();

        Dynamic_Matrix_ solution_delta_matrix(num_equation, num_near_cell);
        for (size_t j = 0; j < num_near_cell; ++j) {
            const auto solution_delta = solutions[near_cell_indexes[j]] - solutions[i];
            for (size_t k = 0; k < num_equation; ++k)
                solution_delta_matrix.at(k, j) = solution_delta[k];
        }
        solution_delta_matrixes.push_back(std::move(solution_delta_matrix));
    }

    return solution_delta_matrixes;
}

template <size_t num_equation, size_t space_dimension>
Vertex_Least_Square<num_equation, space_dimension>::Vertex_Least_Square(const Grid<space_dimension>& grid) {
    SET_TIME_POINT;

    const auto& cell_elements = grid.elements.cell_elements;
    const auto& vnode_index_to_share_cell_indexes = grid.connectivity.vnode_index_to_share_cell_indexes;

    this->num_cell_ = cell_elements.size();
    this->near_cell_indexes_set_.reserve(this->num_cell_);
    this->least_square_matrixes_.reserve(this->num_cell_);

    for (size_t i = 0; i < this->num_cell_; ++i) {
        const auto& element = cell_elements[i];
        const auto& geometry = cell_elements[i].geometry_;

        // near cell indexes - vertex
        auto vnode_indexes = element.vertex_node_indexes();
        std::set<size_t> near_cell_indexes_temp;
        for (const auto vnode_index : vnode_indexes) {
            const auto& share_cell_indexes = vnode_index_to_share_cell_indexes.at(vnode_index);
            near_cell_indexes_temp.insert(share_cell_indexes.begin(), share_cell_indexes.end());
        }
        near_cell_indexes_temp.erase(i);
        std::vector<size_t> near_cell_indexes(near_cell_indexes_temp.begin(), near_cell_indexes_temp.end());

        //least square matrix
        const auto num_neighbor_cell = near_cell_indexes.size();

        const auto this_center = geometry.center_node();

        Dynamic_Matrix_ center_to_center_matrix(space_dimension, num_neighbor_cell);
        for (size_t i = 0; i < num_neighbor_cell; ++i) {
            const auto& neighbor_geometry = cell_elements[near_cell_indexes[i]].geometry_;
            const auto neighbor_center = neighbor_geometry.center_node();
            const auto center_to_center = neighbor_center - this_center;
            for (size_t j = 0; j < space_dimension; ++j)
                center_to_center_matrix.at(j, i) = center_to_center[j];
        }

        const auto& Rc = center_to_center_matrix;
        auto RcT = Rc.transpose();
        auto least_square_matrix = RcT * (Rc * RcT).be_inverse();

        this->near_cell_indexes_set_.push_back(std::move(near_cell_indexes));
        this->least_square_matrixes_.push_back(std::move(least_square_matrix));
    }

    Log::content_ << std::left << std::setw(50) << "@ Vertex Least Sqaure precalculation" << " ----------- " << GET_TIME_DURATION << "s\n\n";
    Log::print();
}



template <size_t num_equation, size_t space_dimension>
Face_Least_Square<num_equation, space_dimension>::Face_Least_Square(const Grid<space_dimension>& grid) {
    SET_TIME_POINT;

    const auto& cell_elements = grid.elements.cell_elements;
    const auto& vnode_index_to_share_cell_indexes = grid.connectivity.vnode_index_to_share_cell_indexes;

    this->num_cell_ = cell_elements.size();
    this->near_cell_indexes_set_.reserve(this->num_cell_);
    this->least_square_matrixes_.reserve(this->num_cell_);

    for (size_t i = 0; i < this->num_cell_; ++i) {
        const auto& element = cell_elements[i];
        const auto& geometry = cell_elements[i].geometry_;

        // near cell indexes - face
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
            this_face_share_cell_indexes.erase(my_index_pos_iter);
            //dynamic_require(this_face_share_cell_indexes.size() == 1, "face share cell should be unique"); // debug
            face_share_cell_indexes.push_back(this_face_share_cell_indexes.front());
        }

        //least square matrix
        const auto num_neighbor_cell = face_share_cell_indexes.size();

        const auto this_center = geometry.center_node();

        Dynamic_Matrix_ center_to_center_matrix(space_dimension, num_neighbor_cell);
        for (size_t i = 0; i < num_neighbor_cell; ++i) {
            const auto& neighbor_geometry = cell_elements[face_share_cell_indexes[i]].geometry_;
            const auto neighbor_center = neighbor_geometry.center_node();
            const auto center_to_center = neighbor_center - this_center;
            for (size_t j = 0; j < space_dimension; ++j)
                center_to_center_matrix.at(j, i) = center_to_center[j];
        }

        const auto& Rc = center_to_center_matrix;
        auto RcT = Rc.transpose();
        auto least_square_matrix = RcT * (Rc * RcT).be_inverse();

        this->near_cell_indexes_set_.push_back(std::move(face_share_cell_indexes));
        this->least_square_matrixes_.push_back(std::move(least_square_matrix));
    }

    Log::content_ << std::left << std::setw(50) << "@ Face Least Sqaure precalculation" << " ----------- " << GET_TIME_DURATION << "s\n\n";
    Log::print();
}