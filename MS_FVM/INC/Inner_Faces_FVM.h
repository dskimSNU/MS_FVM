#pragma once
#include "Grid_Builder.h"
#include "Reconstruction_Method.h"


//FVM이면 공통으로 사용하는 variable
template <size_t space_dimension>
class Inner_Faces_FVM_Base
{
private:
    using Space_Vector_ = Euclidean_Vector<space_dimension>;

protected:
    size_t num_inner_face_ = 0;
    std::vector<Space_Vector_> normals_;
    std::vector<std::pair<size_t, size_t>> oc_nc_index_pairs_;
    std::vector<double> areas_;

public:
    Inner_Faces_FVM_Base(Grid<space_dimension>&& grid);
};


//FVM이고 Constant Reconstruction이면 공통으로 사용하는 variable & Method
template <size_t space_dimension>
class Inner_Faces_FVM_Constant : public Inner_Faces_FVM_Base<space_dimension>
{
private:
    using Space_Vector_ = Euclidean_Vector<space_dimension>;

public:
    Inner_Faces_FVM_Constant(Grid<space_dimension>&& grid) : Inner_Faces_FVM_Base<space_dimension>(std::move(grid)) {};

    template<typename Numerical_Flux_Function, typename Residual, typename Solution>
    void calculate_RHS(std::vector<Residual>& RHS, const std::vector<Solution>& solutions) const;
};



//FVM이고 Linear Reconstruction이면 공통으로 사용하는 variable & Method
template <size_t space_dimension>
class Inner_Faces_FVM_Linear : public Inner_Faces_FVM_Base<space_dimension>
{
private:
    using Space_Vector_ = Euclidean_Vector<space_dimension>;

protected:
    std::vector<std::pair<Space_Vector_, Space_Vector_>> oc_nc_to_face_vector_pairs_;

public:
    Inner_Faces_FVM_Linear(Grid<space_dimension>&& grid);

    template<typename Numerical_Flux_Function, size_t num_equation>
    void calculate_RHS(std::vector<Euclidean_Vector<num_equation>>& RHS, const Linear_Reconstructed_Solution<num_equation, space_dimension>& linear_reconstructed_solution) const;
};


// template definition part
template <size_t space_dimension>
Inner_Faces_FVM_Base<space_dimension>::Inner_Faces_FVM_Base(Grid<space_dimension> && grid) {
    SET_TIME_POINT;

    this->num_inner_face_ = grid.elements.inner_face_elements.size();

    this->areas_.reserve(this->num_inner_face_);
    for (const auto& element : grid.elements.inner_face_elements)
        this->areas_.push_back(element.geometry_.volume());

    this->normals_ = std::move(grid.connectivity.inner_face_normals);
    this->oc_nc_index_pairs_ = std::move(grid.connectivity.inner_face_oc_nc_index_pairs);

    Log::content_ << std::left << std::setw(50) << "@ Inner faces FVM base precalculation" << " ----------- " << GET_TIME_DURATION << "s\n\n";
    Log::print();
}


template <size_t space_dimension>
template<typename Numerical_Flux_Function, typename Residual, typename Solution>
void Inner_Faces_FVM_Constant<space_dimension>::calculate_RHS(std::vector<Residual>& RHS, const std::vector<Solution>& solutions) const {
    const auto numerical_fluxes = Numerical_Flux_Function::calculate(solutions, this->normals_, this->oc_nc_index_pairs_);
    for (size_t i = 0; i < this->num_inner_face_; ++i) {
        const auto [oc_index, nc_index] = this->oc_nc_index_pairs_[i];
        const auto delta_RHS = this->areas_[i] * numerical_fluxes[i];
        RHS[oc_index] -= delta_RHS;
        RHS[nc_index] += delta_RHS;
    }
}


template <size_t space_dimension>
Inner_Faces_FVM_Linear<space_dimension>::Inner_Faces_FVM_Linear(Grid<space_dimension>&& grid) : Inner_Faces_FVM_Base<space_dimension>(std::move(grid)) {
    SET_TIME_POINT;

    this->oc_nc_to_face_vector_pairs_.reserve(this->num_inner_face_);

    const auto& cell_elements = grid.elements.cell_elements;
    for (size_t i = 0; i < this->num_inner_face_; ++i) {
        const auto [oc_index, nc_index] = this->oc_nc_index_pairs_[i];

        const auto& oc_geometry = cell_elements[oc_index].geometry_;
        const auto& nc_geometry = cell_elements[nc_index].geometry_;
        const auto& inner_face_geometry = grid.elements.inner_face_elements[i].geometry_;

        const auto oc_center = oc_geometry.center_node();
        const auto nc_center = nc_geometry.center_node();
        const auto inner_face_center = inner_face_geometry.center_node();

        const auto oc_to_face_vector = inner_face_center - oc_center;
        const auto nc_to_face_vector = inner_face_center - nc_center;

        this->oc_nc_to_face_vector_pairs_.push_back(std::make_pair(oc_to_face_vector, nc_to_face_vector));
    }

    Log::content_ << std::left << std::setw(50) << "@ Inner faces FVM linear precalculation" << " ----------- " << GET_TIME_DURATION << "s\n\n";
    Log::print();
};


template <size_t space_dimension>
template <typename Numerical_Flux_Function, size_t num_equation>
void Inner_Faces_FVM_Linear<space_dimension>::calculate_RHS(std::vector<Euclidean_Vector<num_equation>>& RHS, const Linear_Reconstructed_Solution<num_equation, space_dimension>& linear_reconstructed_solution) const {
    const auto& solutions = linear_reconstructed_solution.solutions;
    const auto& solution_gradients = linear_reconstructed_solution.solution_gradients;

    for (size_t i = 0; i < this->num_inner_face_; ++i) {
        const auto [oc_index, nc_index] = this->oc_nc_index_pairs_[i];
        const auto& oc_solution = solutions[oc_index];
        const auto& nc_solution = solutions[nc_index];

        const auto& oc_solution_gradient = solution_gradients[oc_index];
        const auto& nc_solution_gradient = solution_gradients[nc_index];

        const auto& [oc_to_face_vector, nc_to_face_vector] = this->oc_nc_to_face_vector_pairs_[i];

        const auto oc_side_solution = oc_solution + oc_solution_gradient * oc_to_face_vector;
        const auto nc_side_solution = nc_solution + nc_solution_gradient * nc_to_face_vector;
        const auto& inner_face_normal = this->normals_[i];

        const auto numerical_flux = Numerical_Flux_Function::calculate(oc_side_solution, nc_side_solution, inner_face_normal);
        const auto delta_RHS = this->areas_[i] * numerical_flux;
        RHS[oc_index] -= delta_RHS;
        RHS[nc_index] += delta_RHS;
    }
}