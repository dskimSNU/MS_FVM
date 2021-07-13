#pragma once
#include "Cells.h"
#include "Inner_Faces.h"
#include "TIme_Step_Method.h"

template <typename SDM, typename TSM, typename G, typename NF>
class Semi_Discrete_Equation
{
    static_require(ms::is_spatial_discrete_method<SDM>, "Wrong spatial discrete method");
    static_require(ms::is_governing_equation<G>,        "Wrong governing equation");
    static_require(ms::is_numeirical_flux_function<NF>, "Wrong numerical flux function");
    static_require(ms::is_time_step_method<TSM>,        "Wrong time step method");

    using Spatial_Discrete_Method       = SDM;
    using Governing_Equation            = G;
    using Numerical_Flux_Func           = NF;
    using Time_Step_Method              = TSM;
    using Solution                      = typename Governing_Equation::Solution;
    using Residual                      = typename EuclideanVector<Governing_Equation::num_equation()>;

    static constexpr size_t dimension_          = Governing_Equation::dimension();
    static constexpr double time_step_constant_ = Time_Step_Method::constant();

private:
    Cells<Spatial_Discrete_Method, Governing_Equation> cells_;
    Inner_Faces<Spatial_Discrete_Method, Numerical_Flux_Func> inner_faces_;

public:
    Semi_Discrete_Equation(Grid_Info_Extractor<Spatial_Discrete_Method, dimension_>::Grid_Infos&& grid_info)
        : cells_(std::move(grid_info.cell_informations)), inner_faces_(std::move(grid_info.inner_face_informations)) {};

    double calculate_time_step(const std::vector<Solution>& solutions) const {
        if constexpr (std::is_same_v<Time_Step_Method, CFL<time_step_constant_>>)
            return this->cells_.calculate_time_step(solutions, time_step_constant_);
        else
            return time_step_constant_;
    }

    std::vector<Residual> calculate_RHS(const std::vector<Solution>& solutions) const {
        static const auto num_solution = solutions.size();
        std::vector<Residual> RHS(num_solution);

        if constexpr (std::is_same_v<Spatial_Discrete_Method, FVM>) {
            this->inner_faces_.calculate_RHS(RHS, solutions);
            this->cells_.scale_RHS(RHS);            
        }
        else if constexpr (std::is_same_v<Spatial_Discrete_Method, FVM_Limiter>) {

        }

        return RHS;
    }

    template <typename Initial_Condition>
    void estimate_error(const std::vector<Solution>& computed_solution, const double time)const {
        cells_.estimate_error<Initial_Condition>(computed_solution, time);
    }
};
