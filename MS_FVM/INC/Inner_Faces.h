#pragma once
#include "Grid_Info_Extractor.h"
#include "Numerical_Flux_Function.h"
#include "Spatial_Discrete_Method.h"


template <typename Spatial_Discrete_Method, typename Reconstruction_Method, typename Numerical_Flux_Function>
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
    std::vector<std::pair<size_t, size_t>> cell_container_indexes_o_n_;
    std::vector<double> areas_;

public:
    Inner_Faces_FVM_Base(const Processed_Grid_Data<space_dimension>& processed_grid_data);
};


//Numerical Flux Function에 관계없이 
//FVM이고 MLP계열의 Reconstruction이면 공통으로 사용하는 variable
template <size_t space_dimension>
class Inner_Faces_FVM_MLP : public Inner_Faces_FVM_Base<space_dimension>
{
    using Space_Vector = EuclideanVector<space_dimension>;

protected:
    std::vector<std::pair<Space_Vector, Space_Vector>> cell_to_face_vectors_o_n_;

public:
    Inner_Faces_FVM_MLP(const Processed_Grid_Data<space_dimension>& processed_grid_data);        
};


template<typename Numerical_Flux_Function>
class Inner_Faces<FVM, Constant_Reconstruction, Numerical_Flux_Function> : public Inner_Faces_FVM_Base<Numerical_Flux_Function::Governing_Equation::dimension()>
{
    static_require(ms::is_numeirical_flux_function<Numerical_Flux_Function>, "Wrong numerical flux");

    using Solution_         = typename Numerical_Flux_Function::Solution;
    using Residual_         = typename Numerical_Flux_Function::Numerical_Flux;

    static constexpr size_t space_dimension_ = Numerical_Flux_Function::Governing_Equation::dimension();

public:
    Inner_Faces(const Processed_Grid_Data<space_dimension_>& processed_grid_data)
        : Inner_Faces_FVM_Base<space_dimension_>(processed_grid_data) {};

    void calculate_RHS(std::vector<Residual_>& RHS, const std::vector<Solution_>& solutions) const;
};


template<typename Numerical_Flux_Function, typename Gradient_Method>
class Inner_Faces<FVM, MLP_u1<Gradient_Method>, Numerical_Flux_Function> : public Inner_Faces_FVM_MLP<Numerical_Flux_Function::Governing_Equation::dimension()>
{
    static_require(ms::is_numeirical_flux_function<Numerical_Flux_Function>, "Wrong numerical flux");

    using Solution_     = typename Numerical_Flux_Function::Solution;
    using Residual_     = typename Numerical_Flux_Function::Numerical_Flux;

    static constexpr size_t num_equation_ = Solution_::dimension();
    static constexpr size_t space_dimension_ = Numerical_Flux_Function::Governing_Equation::dimension();

    using Gradient_ = Matrix<num_equation_, space_dimension_>;

public:
    Inner_Faces(const Processed_Grid_Data<space_dimension_>& processed_grid_data)
        : Inner_Faces_FVM_MLP<space_dimension_>(processed_grid_data) {};

    void calculate_RHS(std::vector<Residual_>& RHS, const std::vector<Solution_>& solutions, const std::vector<Gradient_>& solution_gradients_) const;
};


//template definition part
template <size_t space_dimension>
Inner_Faces_FVM_Base<space_dimension>::Inner_Faces_FVM_Base(const Processed_Grid_Data<space_dimension>& processed_grid_data) {
    // periodic boundary can be seen as inner face
    const auto num_periodic_boundary_pair = processed_grid_data.periodic_boundary_owner_side_geometries.size();
    const auto num_inner_face = processed_grid_data.inner_face_geometries.size();
    this->num_inner_face_ = num_periodic_boundary_pair + num_inner_face;

    this->areas_.reserve(this->num_inner_face_);
    this->normals_.reserve(this->num_inner_face_);
    this->cell_container_indexes_o_n_.reserve(this->num_inner_face_);

    const auto& cell_geometries = processed_grid_data.cell_geometries;
    for (size_t i = 0; i < num_periodic_boundary_pair; ++i) {
        const auto& [container_index_o, container_index_n] = processed_grid_data.periodic_boundary_owner_neighbor_container_indexes[i];
        const auto& owner_cell_geometry = cell_geometries.at(container_index_o);
        const auto& owner_side_geometry = processed_grid_data.periodic_boundary_owner_side_geometries[i];

        this->areas_.push_back(owner_side_geometry.volume());
        this->normals_.push_back(owner_side_geometry.normal_vector(owner_cell_geometry.center_node()));
        this->cell_container_indexes_o_n_.push_back({ container_index_o,container_index_n });
    }

    for (size_t i = 0; i < num_inner_face; ++i) {
        const auto& [container_index_o, container_index_n] = processed_grid_data.inner_face_owner_neighbor_container_indexes[i];
        const auto& owner_cell_geometry = cell_geometries.at(container_index_o);
        const auto& inner_face_geometry = processed_grid_data.inner_face_geometries[i];
             
        this->areas_.push_back(inner_face_geometry.volume());
        this->normals_.push_back(inner_face_geometry.normal_vector(owner_cell_geometry.center_node()));
        this->cell_container_indexes_o_n_.push_back({ container_index_o,container_index_n });
    }
}

template <typename Numerical_Flux_Function>
void Inner_Faces<FVM, Constant_Reconstruction, Numerical_Flux_Function>::calculate_RHS(std::vector<Residual_>& RHS, const std::vector<Solution_>& solutions) const {
    const auto numerical_fluxes = Numerical_Flux_Function::calculate(solutions, this->normals_, this->cell_container_indexes_o_n_);
    for (size_t i = 0; i < this->num_inner_face_; ++i) {
        const auto [container_index_o, container_index_n] = this->cell_container_indexes_o_n_[i];
        const auto delta_RHS = this->areas_[i] * numerical_fluxes[i];
        RHS[container_index_o] -= delta_RHS;
        RHS[container_index_n] += delta_RHS;
    }
};

template <size_t space_dimension>
Inner_Faces_FVM_MLP<space_dimension>::Inner_Faces_FVM_MLP(const Processed_Grid_Data<space_dimension>& processed_grid_data)
    :Inner_Faces_FVM_Base<space_dimension>(processed_grid_data) {
    this->cell_to_face_vectors_o_n_.reserve(this->num_inner_face_);

    const auto& cell_geometries = processed_grid_data.cell_geometries;

    //periodic boundary
    const auto num_periodic_boundary_pair = processed_grid_data.periodic_boundary_owner_side_geometries.size();
    for (size_t i = 0; i < num_periodic_boundary_pair; ++i) {
        const auto& [cell_container_index_o, cell_container_index_n] = processed_grid_data.periodic_boundary_owner_neighbor_container_indexes[i];

        const auto& cell_geometry_o = cell_geometries[cell_container_index_o];
        const auto& cell_geometry_n = cell_geometries[cell_container_index_n];
        const auto& owner_side_geometry = processed_grid_data.periodic_boundary_owner_side_geometries[i];

        const auto cell_center_node_o = cell_geometry_o.center_node();
        const auto cell_center_node_n = cell_geometry_n.center_node();
        const auto inner_face_center_node = owner_side_geometry.center_node();

        const auto cell_to_face_vector_o = cell_center_node_o - inner_face_center_node;
        const auto cell_to_face_vector_n = cell_center_node_n - inner_face_center_node;

        this->cell_to_face_vectors_o_n_.push_back({ cell_to_face_vector_o,cell_to_face_vector_n });
     }

    //inner face
    const auto num_inner_face = processed_grid_data.inner_face_geometries.size();
    for (size_t i = 0; i < num_inner_face; ++i) {
        const auto& [cell_container_index_o, cell_container_index_n] = processed_grid_data.inner_face_owner_neighbor_container_indexes[i];

        const auto& cell_geometry_o = cell_geometries[cell_container_index_o];
        const auto& cell_geometry_n = cell_geometries[cell_container_index_n];
        const auto& inner_face_geometry = processed_grid_data.inner_face_geometries[i];

        const auto cell_center_node_o = cell_geometry_o.center_node();
        const auto cell_center_node_n = cell_geometry_n.center_node();
        const auto inner_face_center_node = inner_face_geometry.center_node();

        const auto cell_to_face_vector_o = cell_center_node_o - inner_face_center_node;
        const auto cell_to_face_vector_n = cell_center_node_n - inner_face_center_node;

        this->cell_to_face_vectors_o_n_.push_back({ cell_to_face_vector_o,cell_to_face_vector_n });
    }
};


template<typename Numerical_Flux_Function, typename Gradient_Method>
void Inner_Faces<FVM, MLP_u1<Gradient_Method>, Numerical_Flux_Function>::calculate_RHS(std::vector<Residual_>& RHS, const std::vector<Solution_>& solutions, const std::vector<Gradient_>& solution_gradients) const {
    for (size_t i = 0; i < this->num_inner_face_; ++i) {
        const auto [container_index_o, container_index_n] = this->cell_container_indexes_o_n_[i];
        const auto& solution_o = solutions[container_index_o];
        const auto& solution_n = solutions[container_index_n];

        const auto& solution_gradient_o = solution_gradients[container_index_o];
        const auto& solution_gradient_n = solution_gradients[container_index_n];

        const auto& [cell_to_face_vector_o, cell_to_face_vector_n] = this->cell_to_face_vectors_o_n_[i];

        const auto face_solution_o = solution_o + solution_gradient_o * cell_to_face_vector_o;
        const auto face_solution_n = solution_n + solution_gradient_n * cell_to_face_vector_n;
        const auto& normal = this->normals_[i];

        const auto numerical_flux = Numerical_Flux_Function::calculate(face_solution_o, face_solution_n, normal);
        const auto delta_RHS = this->areas_[i] * numerical_flux;
        RHS[container_index_o] -= delta_RHS;
        RHS[container_index_n] += delta_RHS;
    }
}
