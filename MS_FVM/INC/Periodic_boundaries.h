#pragma once

#pragma once
#include "Grid_Builder.h"
#include "Numerical_Flux_Function.h"
#include "Spatial_Discrete_Method.h"
#include "Reconstruction_Method.h"


template <typename Governing_Equation, typename Spatial_Discrete_Method, typename Reconstruction_Method>
class Periodic_Boundaries;


//Reconstruction Method, Numerical Flux Function에 관계없이
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


////Numerical Flux Function에 관계없이 
////FVM이고 MLP계열의 Reconstruction이면 공통으로 사용하는 variable
//template <size_t space_dimension>
//class Periodic_Boundaries_FVM_MLP_Base : public Periodic_Boundaries_FVM_Base<space_dimension>
//{
//    using Space_Vector = EuclideanVector<space_dimension>;
//
//protected:
//    std::vector<std::pair<Space_Vector, Space_Vector>> oc_nc_to_face_vector_pairs_;
//
//public:
//    Periodic_Boundaries_FVM_MLP_Base(const Grid<space_dimension>& grid);
//};


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


//template<typename Numerical_Flux_Function, typename Gradient_Method>
//class Periodic_Boundaries<FVM, MLP_u1<Gradient_Method>, Numerical_Flux_Function> : public Periodic_Boundaries_FVM_MLP_Base<Numerical_Flux_Function::Governing_Equation::dimension()>
//{
//    static_require(ms::is_numeirical_flux_function<Numerical_Flux_Function>, "Wrong numerical flux");
//
//    using Solution_ = typename Numerical_Flux_Function::Solution;
//    using Residual_ = typename Numerical_Flux_Function::Numerical_Flux;
//
//    static constexpr size_t num_equation_ = Solution_::dimension();
//    static constexpr size_t space_dimension_ = Numerical_Flux_Function::Governing_Equation::dimension();
//
//    using Gradient_ = Matrix<num_equation_, space_dimension_>;
//
//public:
//    Periodic_Boundaries(const Grid<space_dimension_>& grid)
//        : Periodic_Boundaries_FVM_MLP_Base<space_dimension_>(grid) {};
//
//    void calculate_RHS(std::vector<Residual_>& RHS, const std::vector<Solution_>& solutions, const std::vector<Gradient_>& solution_gradients_) const;
//};


//template definition part
template <size_t space_dimension>
Periodic_Boundaries_FVM_Base<space_dimension>::Periodic_Boundaries_FVM_Base(Grid<space_dimension>&& grid) {
    // periodic boundary can be seen as inner face
    this->num_pbdry_pair_ = grid.elements.periodic_boundary_element_pairs.size();

    this->areas_.reserve(this->num_pbdry_pair_);
    for (const auto& [oc_side_element, nc_side_element] : grid.elements.periodic_boundary_element_pairs) {
        this->areas_.push_back(oc_side_element.geometry_.volume());
    }


    this->normals_ = std::move(grid.connectivity.periodic_boundary_normals);
    this->oc_nc_index_pairs_ = std::move(grid.connectivity.periodic_boundary_oc_nc_index_pairs);
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


//template <size_t space_dimension>
//Periodic_Boundaries_FVM_MLP_Base<space_dimension>::Periodic_Boundaries_FVM_MLP_Base(const Grid<space_dimension>& grid)
//    :Periodic_Boundaries_FVM_Base<space_dimension>(grid) {
//    this->oc_nc_to_face_vector_pairs_.reserve(this->num_inner_face_);
//
//    const auto& cell_geometries = grid.cell_geometries;
//
//    //periodic boundary
//    const auto num_periodic_boundary_pair = grid.periodic_boundary_owner_side_geometries.size();
//    for (size_t i = 0; i < num_periodic_boundary_pair; ++i) {
//        const auto& [cell_container_index_o, cell_container_index_n] = grid.periodic_boundary_oc_nc_index_pairs[i];
//
//        const auto& cell_geometry_o = cell_geometries[cell_container_index_o];
//        const auto& cell_geometry_n = cell_geometries[cell_container_index_n];
//        const auto& owner_side_geometry = grid.periodic_boundary_owner_side_geometries[i];
//
//        const auto cell_center_node_o = cell_geometry_o.center_node();
//        const auto cell_center_node_n = cell_geometry_n.center_node();
//        const auto inner_face_center_node = owner_side_geometry.center_node();
//
//        const auto cell_to_face_vector_o = cell_center_node_o - inner_face_center_node;
//        const auto cell_to_face_vector_n = cell_center_node_n - inner_face_center_node;
//
//        this->oc_nc_to_face_vector_pairs_.push_back({ cell_to_face_vector_o,cell_to_face_vector_n });
//    }
//
//    //inner face
//    const auto num_inner_face = grid.inner_face_geometries.size();
//    for (size_t i = 0; i < num_inner_face; ++i) {
//        const auto& [cell_container_index_o, cell_container_index_n] = grid.inner_face_container_indexes_oc_nc[i];
//
//        const auto& cell_geometry_o = cell_geometries[cell_container_index_o];
//        const auto& cell_geometry_n = cell_geometries[cell_container_index_n];
//        const auto& inner_face_geometry = grid.inner_face_geometries[i];
//
//        const auto cell_center_node_o = cell_geometry_o.center_node();
//        const auto cell_center_node_n = cell_geometry_n.center_node();
//        const auto inner_face_center_node = inner_face_geometry.center_node();
//
//        const auto cell_to_face_vector_o = cell_center_node_o - inner_face_center_node;
//        const auto cell_to_face_vector_n = cell_center_node_n - inner_face_center_node;
//
//        this->oc_nc_to_face_vector_pairs_.push_back({ cell_to_face_vector_o,cell_to_face_vector_n });
//    }
//};
//
//
//template<typename Numerical_Flux_Function, typename Gradient_Method>
//void Periodic_Boundaries<FVM, MLP_u1<Gradient_Method>, Numerical_Flux_Function>::calculate_RHS(std::vector<Residual_>& RHS, const std::vector<Solution_>& solutions, const std::vector<Gradient_>& solution_gradients) const {
//    for (size_t i = 0; i < this->num_inner_face_; ++i) {
//        const auto [container_index_o, container_index_n] = this->cell_container_indexes_o_n_[i];
//        const auto& solution_o = solutions[container_index_o];
//        const auto& solution_n = solutions[container_index_n];
//
//        const auto& solution_gradient_o = solution_gradients[container_index_o];
//        const auto& solution_gradient_n = solution_gradients[container_index_n];
//
//        const auto& [cell_to_face_vector_o, cell_to_face_vector_n] = this->cell_to_face_vectors_o_n_[i];
//
//        const auto face_solution_o = solution_o + solution_gradient_o * cell_to_face_vector_o;
//        const auto face_solution_n = solution_n + solution_gradient_n * cell_to_face_vector_n;
//        const auto& normal = this->normals_[i];
//
//        const auto numerical_flux = Numerical_Flux_Function::calculate(face_solution_o, face_solution_n, normal);
//        const auto delta_RHS = this->areas_[i] * numerical_flux;
//        RHS[container_index_o] -= delta_RHS;
//        RHS[container_index_n] += delta_RHS;
//    }
//}
