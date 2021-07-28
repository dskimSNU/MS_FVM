#pragma once
#include "Boundary_Flux_Function.h"
#include "Grid_Builder.h"
#include "Reconstruction_Method.h"


//FVM이면 공통으로 사용하는 variable
template <typename Governing_Equation>
class Boundaries_FVM_Base
{
private:
    static constexpr size_t space_dimension_ = Governing_Equation::space_dimension();

    using Space_Vector_ = Governing_Equation::Space_Vector_;

protected:
    size_t num_boundaries_ = 0;
    std::vector<Space_Vector_> normals_;
    std::vector<size_t> oc_indexes_;
    std::vector<double> areas_;    
    std::vector<std::unique_ptr<Boundary_Flux_Function<Governing_Equation>>> boundary_flux_functions_;

public:
    Boundaries_FVM_Base(Grid<space_dimension_>&& grid);
};


//FVM이고 Constant Reconstruction이면 공통으로 사용하는 variable & Method
template <typename Governing_Equation>
class Boundaries_FVM_Constant : public Boundaries_FVM_Base<Governing_Equation>
{
private:
    static constexpr size_t space_dimension_ = Governing_Equation::space_dimension();
    static constexpr size_t num_equation_ = Governing_Equation::num_equation();

    using Solution_     = typename Governing_Equation::Solution_;
    using Boundary_Flux_     = EuclideanVector<num_equation_>;

public:
    Boundaries_FVM_Constant(Grid<space_dimension_>&& grid) : Boundaries_FVM_Base<Governing_Equation>(std::move(grid)) {};

    void calculate_RHS(std::vector<Boundary_Flux_>& RHS, const std::vector<Solution_>& solutions) const;
};


//FVM이고 Constant Reconstruction이면 공통으로 사용하는 variable & Method
template <typename Governing_Equation>
class Boundaries_FVM_Linear : public Boundaries_FVM_Base<Governing_Equation>
{
private:
    static constexpr size_t space_dimension_ = Governing_Equation::space_dimension();
    static constexpr size_t num_equation_ = Governing_Equation::num_equation();

    using Space_Vector_ = EuclideanVector <space_dimension_>;
    using Solution_     = typename Governing_Equation::Solution_;
    using Boundary_Flux_     = EuclideanVector<num_equation_>;

private:
    std::vector<Space_Vector_> oc_to_boundary_vectors_;

public:
    Boundaries_FVM_Linear(Grid<space_dimension_>&& grid);

    void calculate_RHS(std::vector<Boundary_Flux_>& RHS, const Linear_Reconstructed_Solution<num_equation_, space_dimension_>& linear_reconstructed_solution) const;
};


//template definition
template <typename Governing_Equation>
Boundaries_FVM_Base<Governing_Equation>::Boundaries_FVM_Base(Grid<space_dimension_>&& grid) {
    SET_TIME_POINT;

    this->num_boundaries_ = grid.elements.boundary_elements.size();

    this->areas_.reserve(this->num_boundaries_);
    this->boundary_flux_functions_.reserve(this->num_boundaries_);
    for (const auto& element : grid.elements.boundary_elements) {
        this->areas_.push_back(element.geometry_.volume());
        this->boundary_flux_functions_.push_back(Boundary_Flux_Function_Factory<Governing_Equation>::make(element.type()));
    }

    this->normals_ = std::move(grid.connectivity.boundary_normals);
    this->oc_indexes_ = std::move(grid.connectivity.boundary_oc_indexes);

    Log::content_ << std::left << std::setw(50) << "@ Boundaries FVM base precalculation" << " ----------- " << GET_TIME_DURATION << "s\n\n";
    Log::print();
}

template <typename Governing_Equation>
void Boundaries_FVM_Constant<Governing_Equation>::calculate_RHS(std::vector<Boundary_Flux_>& RHS, const std::vector<Solution_>& solutions) const {
    for (size_t i = 0; i < this->num_boundaries_; ++i) {        
        const auto& boundary_flux_function = this->boundary_flux_functions_.at(i);
        
        const auto oc_index = this->oc_indexes_[i];
        const auto& normal = this->normals_[i];

        const auto boundary_flux = boundary_flux_function->calculate(solutions[oc_index], normal);
        const auto delta_RHS = this->areas_[i] * boundary_flux;

        RHS[oc_index] -= delta_RHS;
    }
}

template <typename Governing_Equation>
Boundaries_FVM_Linear<Governing_Equation>::Boundaries_FVM_Linear(Grid<space_dimension_>&& grid) : Boundaries_FVM_Base<Governing_Equation>(std::move(grid)) {
    SET_TIME_POINT;

    this->oc_to_boundary_vectors_.reserve(this->num_boundaries_);

    const auto& cell_elements = grid.elements.cell_elements;
    const auto& boundary_elements = grid.elements.boundary_elements;
    for (size_t i = 0; i < this->num_boundaries_; ++i) {
        const auto oc_index = this->oc_indexes_[i];

        const auto& oc_geometry = cell_elements[oc_index].geometry_;
        const auto& boundary_geometry = boundary_elements[i].geometry_;

        const auto oc_center = oc_geometry.center_node();
        const auto boundary_center = boundary_geometry.center_node();

        const auto oc_to_face_vector = boundary_center - oc_center;

        this->oc_to_boundary_vectors_.push_back(oc_to_face_vector);
    }

    Log::content_ << std::left << std::setw(50) << "@ Boundaries FVM linear precalculation" << " ----------- " << GET_TIME_DURATION << "s\n\n";
    Log::print();
};

template <typename Governing_Equation>
void Boundaries_FVM_Linear<Governing_Equation>::calculate_RHS(std::vector<Boundary_Flux_>& RHS, const Linear_Reconstructed_Solution<num_equation_, space_dimension_>& linear_reconstructed_solution) const {
    const auto& solutions = linear_reconstructed_solution.solutions;
    const auto& solution_gradients = linear_reconstructed_solution.solution_gradients;

    for (size_t i = 0; i < this->num_boundaries_; ++i) {
        const auto& boundary_flux_function = this->boundary_flux_functions_.at(i);
        const auto& normal = this->normals_.at(i);

        const auto oc_index = this->oc_indexes_.at(i);

        const auto& oc_solution = solutions[oc_index];
        const auto& oc_solution_gradient = solution_gradients[oc_index];
        const auto& oc_to_face_vector = this->oc_to_boundary_vectors_[i];

        const auto oc_side_solution = oc_solution + oc_solution_gradient * oc_to_face_vector;

        const auto boundary_flux = boundary_flux_function->calculate(oc_side_solution, normal);

        const auto delta_RHS = this->areas_[i] * boundary_flux;

        RHS[oc_index] -= delta_RHS;
    }
}