#pragma once

#include "Cells.h"
#include "Inner_Faces.h"

class Semi_Discrete_Equation {};

template<typename G, typename NF>
class FVM_Semi_Discrete_Equation : Semi_Discrete_Equation
{
    static_require(ms::is_Gov_Eq<G>, "Wrong governing equation");
    static_require(ms::is_Numerical_Flux_Func<NF>, "Wrong numerical flux function");

public:
    using Gov_Eq = G;
    using Numerical_Flux_Func = NF;
    using Solution = Gov_Eq::Solution;
    using Residual = EuclideanVector<Gov_Eq::num_equation_>;

public:
    FVM_Semi_Discrete_Equation(Grid_Information&& grid_info)
        : cells_(std::move(grid_info.cell_grid_information)), inner_faces_(std::move(grid_info.inner_face_grid_information)) {};

    double calculate_time_step(const std::vector<Solution>& solutions) const;
    std::vector<Residual> calculate_RHS(const std::vector<Solution>& solutions) const;

private:
    Cells<Gov_Eq> cells_;
    Inner_Faces<Numerical_Flux_Func> inner_faces_;
};

namespace ms{
    template <typename T>
    inline constexpr bool is_Semi_Discrete_Eq = std::is_base_of_v<Semi_Discrete_Equation, T>;
}

//template definition
template<typename G, typename NF>
std::vector<typename FVM_Semi_Discrete_Equation<G,NF>::Residual> FVM_Semi_Discrete_Equation<G,NF>::calculate_RHS(const std::vector<Solution>& solutions) const {
    static const auto num_solution = solutions.size();
    std::vector<Residual> RHS(num_solution);
    this->inner_faces_.calculate_RHS(RHS, solutions);
    this->cells_.scale_RHS(RHS);
    return RHS;
}

template<typename G, typename NF>
double FVM_Semi_Discrete_Equation<G, NF>::calculate_time_step(const std::vector<Solution>& solutions) const {
    return this->cells_.calculate_time_step(solutions);
}