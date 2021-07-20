#pragma once

#include "Grid_Builder.h"
#include "Spatial_Discrete_Method.h"
#include "Reconstruction_Method.h"


//FVM이면 공통으로 사용하는 variable
template <size_t space_dimension>
class Periodic_Boundaries_FVM_Base
{
    using Space_Vector = EuclideanVector<space_dimension>;

protected:
    size_t num_pbdry_pair_ = 0;
    std::vector<Space_Vector> normals_;
    std::vector<std::pair<size_t, size_t>> oc_nc_index_pairs_;
    std::vector<double> areas_;

public:
    Periodic_Boundaries_FVM_Base(Grid<space_dimension>&& grid);
};


//FVM이고 Linear Reconstruction이면 공통으로 사용하는 variable
template <size_t space_dimension>
class Periodic_Boundaries_FVM_Linear_Base : public Periodic_Boundaries_FVM_Base<space_dimension>
{
private:
    using Space_Vector = EuclideanVector<space_dimension>;

protected:
    std::vector<std::pair<Space_Vector, Space_Vector>> oc_nc_to_oc_nc_side_face_vector_pairs_;

public:
    Periodic_Boundaries_FVM_Linear_Base(Grid<space_dimension>&& grid);

    template<typename Numerical_Flux_Function, typename Residual, typename Solution, typename Gradient>
    void calculate_RHS_with_gradient(std::vector<Residual>& RHS, const std::vector<Solution>& solutions, const std::vector<Gradient>& solution_gradients) const;
};
////Numerical Flux Function에 관계없이 
////FVM이고 MLP계열의 Reconstruction이면 공통으로 사용하는 variable
//template <size_t space_dimension>
//class Periodic_Boundaries_FVM_MLP_Base : public Periodic_Boundaries_FVM_Base<space_dimension>
//{
//    using Space_Vector = EuclideanVector<space_dimension>;
//
//protected:
//    std::vector<std::pair<Space_Vector, Space_Vector>> oc_nc_to_oc_nc_side_face_vector_pairs_;
//
//public:
//    Periodic_Boundaries_FVM_MLP_Base(Grid<space_dimension>&& grid);
//};


template <typename Governing_Equation, typename Spatial_Discrete_Method, typename Reconstruction_Method>
class Periodic_Boundaries;


template<typename Governing_Equation>
class Periodic_Boundaries<Governing_Equation, FVM, Constant_Reconstruction> : public Periodic_Boundaries_FVM_Base<Governing_Equation::space_dimension()>
{
private:
    using Solution_ = typename Governing_Equation::Solution_;
    using Residual_ = EuclideanVector<Governing_Equation::num_equation()>;

    static constexpr size_t space_dimension_ = Governing_Equation::space_dimension();

public:
    Periodic_Boundaries(Grid<space_dimension_>&& grid) : Periodic_Boundaries_FVM_Base<space_dimension_>(std::move(grid)) {};

    template <typename Numerical_Flux_Function>
    void calculate_RHS(std::vector<Residual_>& RHS, const std::vector<Solution_>& solutions) const;
};

template<typename Governing_Equation, typename Gradient_Method>
class Periodic_Boundaries<Governing_Equation, FVM, MLP_u1<Gradient_Method>> : public Periodic_Boundaries_FVM_Linear_Base<Governing_Equation::space_dimension()>
{
private:
    static constexpr size_t space_dimension_ = Governing_Equation::space_dimension();

public:
    Periodic_Boundaries(Grid<space_dimension_>&& grid) : Periodic_Boundaries_FVM_Linear_Base<space_dimension_>(std::move(grid)) {};
};


//template<typename Governing_Equation, typename Gradient_Method>
//class Periodic_Boundaries<Governing_Equation, FVM, MLP_u1<Gradient_Method>> : public Periodic_Boundaries_FVM_MLP_Base<Governing_Equation::space_dimension()>
//{
//private:
//    static constexpr size_t space_dimension_ = Governing_Equation::space_dimension();
//    static constexpr size_t num_equation_ = Governing_Equation::num_equation();
//
//    using Solution_ = typename Governing_Equation::Solution_;
//    using Residual_ = EuclideanVector<num_equation_>;
//    using Gradient_ = Matrix<num_equation_, space_dimension_>;
//
//public:
//    Periodic_Boundaries(Grid<space_dimension_>&& grid) : Periodic_Boundaries_FVM_MLP_Base<space_dimension_>(std::move(grid)) {};
//
//    template <typename Numerical_Flux_Function>
//    void calculate_RHS(std::vector<Residual_>& RHS, const std::vector<Solution_>& solutions, const std::vector<Gradient_>& solution_gradients_) const;
//};


//template definition part
template <size_t space_dimension>
Periodic_Boundaries_FVM_Base<space_dimension>::Periodic_Boundaries_FVM_Base(Grid<space_dimension>&& grid) {
    SET_TIME_POINT;
    
    this->num_pbdry_pair_ = grid.elements.periodic_boundary_element_pairs.size();

    this->areas_.reserve(this->num_pbdry_pair_);
    for (const auto& [oc_side_element, nc_side_element] : grid.elements.periodic_boundary_element_pairs) {
        this->areas_.push_back(oc_side_element.geometry_.volume());
    }

    this->normals_ = std::move(grid.connectivity.periodic_boundary_normals);
    this->oc_nc_index_pairs_ = std::move(grid.connectivity.periodic_boundary_oc_nc_index_pairs);

    Log::content_ << std::left << std::setw(50) << "@ Construct Periodic boundaries FVM Base" << " ----------- " << GET_TIME_DURATION << "s\n\n";
    Log::print();
}

template <typename Governing_Equation>
template <typename Numerical_Flux_Function>
void Periodic_Boundaries<Governing_Equation, FVM, Constant_Reconstruction>::calculate_RHS(std::vector<Residual_>& RHS, const std::vector<Solution_>& solutions) const {
    const auto numerical_fluxes = Numerical_Flux_Function::calculate(solutions, this->normals_, this->oc_nc_index_pairs_);
    for (size_t i = 0; i < this->num_pbdry_pair_; ++i) {
        const auto [oc_index, nc_index] = this->oc_nc_index_pairs_[i];
        const auto delta_RHS = this->areas_[i] * numerical_fluxes[i];
        RHS[oc_index] -= delta_RHS;
        RHS[nc_index] += delta_RHS;
    }
};

template <size_t space_dimension>
Periodic_Boundaries_FVM_Linear_Base<space_dimension>::Periodic_Boundaries_FVM_Linear_Base(Grid<space_dimension>&& grid) :Periodic_Boundaries_FVM_Base<space_dimension>(std::move(grid)) {
    SET_TIME_POINT;

    this->oc_nc_to_oc_nc_side_face_vector_pairs_.reserve(this->num_pbdry_pair_);

    const auto& cell_elements = grid.elements.cell_elements;
    const auto& pbdry_element_pairs = grid.elements.periodic_boundary_element_pairs;
    for (size_t i = 0; i < this->num_pbdry_pair_; ++i) {
        const auto& [oc_index, nc_index] = this->oc_nc_index_pairs_[i];
        const auto& oc_geometry = cell_elements[oc_index].geometry_;
        const auto& nc_geometry = cell_elements[nc_index].geometry_;

        const auto& [oc_side_element, nc_side_element] = pbdry_element_pairs[i];
        const auto& oc_side_geometry = oc_side_element.geometry_;
        const auto& nc_side_geometry = nc_side_element.geometry_;

        const auto oc_center = oc_geometry.center_node();
        const auto nc_center = nc_geometry.center_node();
        const auto oc_side_center = oc_side_geometry.center_node();
        const auto nc_side_center = nc_side_geometry.center_node();

        const auto oc_to_oc_side_face_vector = oc_center - oc_side_center;
        const auto nc_to_nc_side_face_vector = nc_center - nc_side_center;

        this->oc_nc_to_oc_nc_side_face_vector_pairs_.push_back(std::make_pair(oc_to_oc_side_face_vector, nc_to_nc_side_face_vector));
    }

    Log::content_ << std::left << std::setw(50) << "@ Construct Periodic boundaries FVM MLP Base" << " ----------- " << GET_TIME_DURATION << "s\n\n";
    Log::print();
}


template <size_t space_dimension>
template<typename Numerical_Flux_Function, typename Residual, typename Solution, typename Gradient>
void Periodic_Boundaries_FVM_Linear_Base<space_dimension>::calculate_RHS_with_gradient(std::vector<Residual>& RHS, const std::vector<Solution>& solutions, const std::vector<Gradient>& solution_gradients) const {
    for (size_t i = 0; i < this->num_pbdry_pair_; ++i) {
        const auto [oc_index, nc_index] = this->oc_nc_index_pairs_[i];
        const auto& oc_solution = solutions[oc_index];
        const auto& nc_solution = solutions[nc_index];

        const auto& oc_solution_gradient = solution_gradients[oc_index];
        const auto& nc_solution_gradient = solution_gradients[nc_index];

        const auto& [oc_to_oc_side_face_vector, nc_to_nc_side_face_vector] = this->oc_nc_to_oc_nc_side_face_vector_pairs_[i];

        const auto oc_side_solution = oc_solution + oc_solution_gradient * oc_to_oc_side_face_vector;
        const auto nc_side_solution = nc_solution + nc_solution_gradient * nc_to_nc_side_face_vector;
        const auto& pbdry_normal = this->normals_[i];

        const auto numerical_flux = Numerical_Flux_Function::calculate(oc_side_solution, nc_side_solution, pbdry_normal);
        const auto delta_RHS = this->areas_[i] * numerical_flux;
        RHS[oc_index] -= delta_RHS;
        RHS[nc_index] += delta_RHS;
    }
}


//template <size_t space_dimension>
//Periodic_Boundaries_FVM_MLP_Base<space_dimension>::Periodic_Boundaries_FVM_MLP_Base(Grid<space_dimension>&& grid) :Periodic_Boundaries_FVM_Base<space_dimension>(std::move(grid)) {
//    SET_TIME_POINT;
//    
//    this->oc_nc_to_oc_nc_side_face_vector_pairs_.reserve(this->num_pbdry_pair_);
//
//    const auto& cell_elements = grid.elements.cell_elements;
//    const auto& pbdry_element_pairs = grid.elements.periodic_boundary_element_pairs;
//    for (size_t i = 0; i < this->num_pbdry_pair_; ++i) {
//        const auto& [oc_index, nc_index] = this->oc_nc_index_pairs_[i];
//        const auto& oc_geometry = cell_elements[oc_index].geometry_;
//        const auto& nc_geometry = cell_elements[nc_index].geometry_;
//
//        const auto& [oc_side_element, nc_side_element] = pbdry_element_pairs[i];
//        const auto& oc_side_geometry = oc_side_element.geometry_;
//        const auto& nc_side_geometry = nc_side_element.geometry_;
//
//        const auto oc_center = oc_geometry.center_node();
//        const auto nc_center = nc_geometry.center_node();
//        const auto oc_side_center = oc_side_geometry.center_node();
//        const auto nc_side_center = nc_side_geometry.center_node();
//
//        const auto oc_to_oc_side_face_vector = oc_center - oc_side_center;
//        const auto nc_to_nc_side_face_vector = nc_center - nc_side_center;
//
//        this->oc_nc_to_oc_nc_side_face_vector_pairs_.push_back(std::make_pair(oc_to_oc_side_face_vector, nc_to_nc_side_face_vector));
//    }
//
//    Log::content_ << std::left << std::setw(50) << "@ Construct Periodic boundaries FVM MLP Base" << " ----------- " << GET_TIME_DURATION << "s\n\n";
//    Log::print();
//};


//template<typename Governing_Equation, typename Gradient_Method>
//template<typename Numerical_Flux_Function>
//void Periodic_Boundaries<Governing_Equation, FVM, MLP_u1<Gradient_Method>>::calculate_RHS(std::vector<Residual_>& RHS, const std::vector<Solution_>& solutions, const std::vector<Gradient_>& solution_gradients) const {
//    for (size_t i = 0; i < this->num_pbdry_pair_; ++i) {
//        const auto [oc_index, nc_index] = this->oc_nc_index_pairs_[i];
//        const auto& oc_solution = solutions[oc_index];
//        const auto& nc_solution = solutions[nc_index];
//
//        const auto& oc_solution_gradient = solution_gradients[oc_index];
//        const auto& nc_solution_gradient = solution_gradients[nc_index];
//
//        const auto& [oc_to_oc_side_face_vector, nc_to_nc_side_face_vector] = this->oc_nc_to_oc_nc_side_face_vector_pairs_[i];
//
//        const auto oc_side_solution = oc_solution + oc_solution_gradient * oc_to_oc_side_face_vector;
//        const auto nc_side_solution = nc_solution + nc_solution_gradient * nc_to_nc_side_face_vector;
//        const auto& pbdry_normal = this->normals_[i];
//
//        const auto numerical_flux = Numerical_Flux_Function::calculate(oc_side_solution, nc_side_solution, pbdry_normal);
//        const auto delta_RHS = this->areas_[i] * numerical_flux;
//        RHS[oc_index] -= delta_RHS;
//        RHS[nc_index] += delta_RHS;
//    }
//}
