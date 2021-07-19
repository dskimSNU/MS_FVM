#pragma once
#include "Cells.h"
#include "Inner_Faces.h"
#include "Periodic_boundaries.h"
#include "Numerical_Flux_Function.h"

template <typename Governing_Equation, typename Spatial_Discrete_Method, typename Reconstruction_Method, typename Numerical_Flux_Function>
class Semi_Discrete_Equation
{
    static_require(ms::is_governing_equation<Governing_Equation>,               "Wrong governing equation");
    static_require(ms::is_spatial_discrete_method<Spatial_Discrete_Method>,     "Wrong spatial discrete method");
    static_require(ms::is_reconsturction_method<Reconstruction_Method>,         "Wrong reconstruction method");
    static_require(ms::is_numeirical_flux_function<Numerical_Flux_Function>,    "Wrong numerical flux function");

    using Cells_                = Cells<Governing_Equation, Spatial_Discrete_Method, Reconstruction_Method>;
    using Periodic_Boundaries_  = Periodic_Boundaries<Governing_Equation, Spatial_Discrete_Method, Reconstruction_Method>;
    using Inner_Faces_          = Inner_Faces<Governing_Equation, Spatial_Discrete_Method, Reconstruction_Method>;
    using Solution_             = typename Governing_Equation::Solution_;
    using Residual_             = EuclideanVector<Governing_Equation::num_equation()>;
    
    static constexpr size_t space_dimension_ = Governing_Equation::space_dimension();

private:
    Cells_ cells_;
    Periodic_Boundaries_ periodic_boundaries_;
    Inner_Faces_ inner_faces_;

public:
    Semi_Discrete_Equation(Grid<space_dimension_>&& grid)
        : cells_(std::move(grid)), periodic_boundaries_(std::move(grid)), inner_faces_(std::move(grid)) {

        Log::content_ << "================================================================================\n";
        Log::content_ << "\t\t\t Total ellapsed time: " << GET_TIME_DURATION << "s\n";
        Log::content_ << "================================================================================\n\n";
        Log::print();
    };


    template <typename Time_Step_Method>
    double calculate_time_step(const std::vector<Solution_>& solutions) const {
        static constexpr double time_step_constant_ = Time_Step_Method::constant();
        if constexpr (std::is_same_v<Time_Step_Method, CFL<time_step_constant_>>) {
            const auto projected_maximum_lambdas = Governing_Equation::coordinate_projected_maximum_lambdas(solutions);
            return this->cells_.calculate_time_step(projected_maximum_lambdas, time_step_constant_);
        }
        else
            return time_step_constant_;
    }

    std::vector<Residual_> calculate_RHS(const std::vector<Solution_>& solutions) const {
        static const auto num_solution = solutions.size();
        std::vector<Residual_> RHS(num_solution);

        if constexpr (std::is_same_v<Reconstruction_Method, Constant_Reconstruction>) {
            this->periodic_boundaries_.calculate_RHS<Numerical_Flux_Function>(RHS, solutions);
            this->inner_faces_.calculate_RHS<Numerical_Flux_Function>(RHS, solutions);
            this->cells_.scale_RHS(RHS);            
        }
        else {
            const auto solution_gradients = this->cells_.calculate_gradient(solutions);
            this->periodic_boundaries_.calculate_RHS<Numerical_Flux_Function>(RHS, solutions, solution_gradients);
            this->inner_faces_.calculate_RHS<Numerical_Flux_Function>(RHS, solutions, solution_gradients);
            this->cells_.scale_RHS(RHS);
        }

        return RHS;
    }

    template <typename Initial_Condition>
    std::vector<Solution_> calculate_initial_solutions(void)const {
        return cells_.calculate_initial_solutions<Initial_Condition>();
    }

    template <typename Initial_Condition>
    void estimate_error(const std::vector<Solution_>& computed_solution, const double time)const {
        cells_.estimate_error<Initial_Condition, Governing_Equation>(computed_solution, time);
    }

};