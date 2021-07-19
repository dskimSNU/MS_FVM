#pragma once
#include "Grid_Builder.h"
#include "Spatial_Discrete_Method.h"
#include "Reconstruction_Method.h"

template <typename Governing_Equation, typename Spatial_Discrete_Method, typename Reconstruction_Method>
class Inner_Faces;


//Reconstruction Method, Numerical Flux Function에 관계없이
//FVM이면 공통으로 사용하는 variable
template <size_t space_dimension>
class Inner_Faces_FVM_Base
{
    using Space_Vector = EuclideanVector<space_dimension>;

protected:
    size_t num_inner_face_ = 0;
    std::vector<Space_Vector> normals_;
    std::vector<std::pair<size_t, size_t>> oc_nc_index_pairs_;
    std::vector<double> areas_;

public:
    Inner_Faces_FVM_Base(Grid<space_dimension>&& grid);
};


//Numerical Flux Function에 관계없이 
//FVM이고 MLP계열의 Reconstruction이면 공통으로 사용하는 variable
template <size_t space_dimension>
class Inner_Faces_FVM_MLP : public Inner_Faces_FVM_Base<space_dimension>
{
    using Space_Vector = EuclideanVector<space_dimension>;

protected:
    std::vector<std::pair<Space_Vector, Space_Vector>> oc_nc_to_face_vector_pairs_;

public:
    Inner_Faces_FVM_MLP(Grid<space_dimension>&& grid);        
};


template<typename Governing_Equation>
class Inner_Faces<Governing_Equation, FVM, Constant_Reconstruction> : public Inner_Faces_FVM_Base<Governing_Equation::space_dimension()>
{
private:
    using Solution_ = typename Governing_Equation::Solution_;
    using Residual_ = EuclideanVector<Governing_Equation::num_equation()>;

    static constexpr size_t space_dimension_ = Governing_Equation::space_dimension();

public:
    Inner_Faces(Grid<space_dimension_>&& grid) : Inner_Faces_FVM_Base<space_dimension_>(std::move(grid)) {};

    template<typename Numerical_Flux_Function>
    void calculate_RHS(std::vector<Residual_>& RHS, const std::vector<Solution_>& solutions) const;
};


template<typename Governing_Equation, typename Gradient_Method>
class Inner_Faces<Governing_Equation, FVM, MLP_u1<Gradient_Method>> : public Inner_Faces_FVM_MLP<Governing_Equation::space_dimension()>
{
private:
    static constexpr size_t space_dimension_ = Governing_Equation::space_dimension();
    static constexpr size_t num_equation_ = Governing_Equation::num_equation();

    using Solution_ = typename Governing_Equation::Solution_;
    using Residual_ = EuclideanVector<num_equation_>;
    using Gradient_ = Matrix<num_equation_, space_dimension_>;

public:
    Inner_Faces(Grid<space_dimension_>&& grid) : Inner_Faces_FVM_MLP<space_dimension_>(std::move(grid)) {};

    template<typename Numerical_Flux_Function>
    void calculate_RHS(std::vector<Residual_>& RHS, const std::vector<Solution_>& solutions, const std::vector<Gradient_>& solution_gradients_) const;
};


//template definition part
template <size_t space_dimension>
Inner_Faces_FVM_Base<space_dimension>::Inner_Faces_FVM_Base(Grid<space_dimension>&& grid) {
    SET_TIME_POINT;

    this->num_inner_face_ = grid.elements.inner_face_elements.size();

    this->areas_.reserve(this->num_inner_face_);
    for (const auto& element : grid.elements.inner_face_elements) 
        this->areas_.push_back(element.geometry_.volume());

    this->normals_ = std::move(grid.connectivity.inner_face_normals);
    this->oc_nc_index_pairs_ = std::move(grid.connectivity.inner_face_oc_nc_index_pairs);

    Log::content_ << std::left << std::setw(50) << "@ Construct Inner faces FVM Base" << " ----------- " << GET_TIME_DURATION << "s\n\n";
    Log::print();
}

template <typename Governing_Equation>
template <typename Numerical_Flux_Function>
void Inner_Faces<Governing_Equation, FVM, Constant_Reconstruction>::calculate_RHS(std::vector<Residual_>& RHS, const std::vector<Solution_>& solutions) const {
    const auto numerical_fluxes = Numerical_Flux_Function::calculate(solutions, this->normals_, this->oc_nc_index_pairs_);
    for (size_t i = 0; i < this->num_inner_face_; ++i) {
        const auto [oc_index, nc_index] = this->oc_nc_index_pairs_[i];
        const auto delta_RHS = this->areas_[i] * numerical_fluxes[i];
        RHS[oc_index] -= delta_RHS;
        RHS[nc_index] += delta_RHS;
    }
};

template <size_t space_dimension>
Inner_Faces_FVM_MLP<space_dimension>::Inner_Faces_FVM_MLP(Grid<space_dimension>&& grid) : Inner_Faces_FVM_Base<space_dimension>(std::move(grid)) {
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

        const auto oc_to_face_vector = oc_center - inner_face_center;
        const auto nc_to_face_vector = nc_center - inner_face_center;

        this->oc_nc_to_face_vector_pairs_.push_back(std::make_pair(oc_to_face_vector, nc_to_face_vector));
    }

    Log::content_ << std::left << std::setw(50) << "@ Construct Inner faces FVM MLP Base" << " ----------- " << GET_TIME_DURATION << "s\n\n";
    Log::print();
};


template<typename Governing_Equation, typename Gradient_Method>
template<typename Numerical_Flux_Function>
void Inner_Faces<Governing_Equation, FVM, MLP_u1<Gradient_Method>>::calculate_RHS(std::vector<Residual_>& RHS, const std::vector<Solution_>& solutions, const std::vector<Gradient_>& solution_gradients) const {
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
