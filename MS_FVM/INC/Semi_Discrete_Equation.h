#pragma once

#include "Cells.h"
#include "Inner_Faces.h"
#include "Spatial_Discrete_Method.h"

template <typename SDM, typename G, typename NF>
class Semi_Discrete_Equation
{
    static_require(ms::is_spatial_discrete_method<SDM>, "Wrong spatial discrete method");
    static_require(ms::is_Gov_Eq<G>, "Wrong governing equation");
    static_require(ms::is_Numerical_Flux_Func<NF>, "Wrong numerical flux function");

public:
    using Spatial_Discrete_Method = SDM;
    using Gov_Eq = G;
    using Numerical_Flux_Func = NF;
    using Solution = Gov_Eq::Solution;
    using Residual = EuclideanVector<Gov_Eq::num_equation_>;

public:
    Semi_Discrete_Equation(Grid_Data_to_Info<Gov_Eq::dimension_>::Grid_Info&& grid_info)
        : cells_(std::move(grid_info.cell_grid_information)), inner_faces_(std::move(grid_info.inner_face_grid_information)) {};

    double calculate_time_step(const std::vector<Solution>& solutions) const {
        return this->cells_.calculate_time_step(solutions);
    }

    std::vector<Residual> calculate_RHS(const std::vector<Solution>& solutions) const {
        static const auto num_solution = solutions.size();
        std::vector<Residual> RHS(num_solution);
        this->inner_faces_.calculate_RHS(RHS, solutions);
        this->cells_.scale_RHS(RHS);
        return RHS;
    }

private:
    Cells<Spatial_Discrete_Method, Gov_Eq> cells_;
    Inner_Faces<Spatial_Discrete_Method, Numerical_Flux_Func> inner_faces_;
};
