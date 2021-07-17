#pragma once
#include "Cells.h"
#include "Inner_Faces.h"
#include "Time_Step_Method.h"

template <typename Governing_Equation, typename Spatial_Discrete_Method, typename Reconstruction_Method, typename Numerical_Flux_Function>
class Semi_Discrete_Equation
{
    static_require(ms::is_governing_equation<Governing_Equation>,               "Wrong governing equation");
    static_require(ms::is_spatial_discrete_method<Spatial_Discrete_Method>,     "Wrong spatial discrete method");
    static_require(ms::is_reconsturction_method<Reconstruction_Method>,         "Wrong reconstruction method");
    static_require(ms::is_numeirical_flux_function<Numerical_Flux_Function>,    "Wrong numerical flux function");

    using Cells_        = typename Cells<Governing_Equation, Spatial_Discrete_Method, Reconstruction_Method>;
    using Inner_Faces_  = typename Inner_Faces<Spatial_Discrete_Method, Reconstruction_Method, Numerical_Flux_Function>;
    using Solution_     = typename Governing_Equation::Solution;
    using Residual_     = typename EuclideanVector<Governing_Equation::num_equation()>;
    
    static constexpr size_t space_dimension_ = Governing_Equation::dimension();

private:
    Cells_ cells_;
    Inner_Faces_ inner_faces_;

public:
    Semi_Discrete_Equation(const Grid<space_dimension_>& grid)
        : cells_(grid), inner_faces_(grid) {};

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
            this->inner_faces_.calculate_RHS(RHS, solutions);
            this->cells_.scale_RHS(RHS);            
        }
        else {
            const auto solution_gradients = this->cells_.calculate_gradient(solutions);
            this->inner_faces_.calculate_RHS(RHS, solutions, solution_gradients);
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