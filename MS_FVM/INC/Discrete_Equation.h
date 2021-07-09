#pragma once
#include "Cells.h"

template<typename NF>
class Semi_Discrete_Equation
{
    using Numerical_Flux_Function = NF;
    using Solution = NF::Solution;
    using Residual = EuclideanVector<NF::Governing_Equation::num_eqation_>;
public:
    Semi_Discrete_Equation(const Cells& cells, Inner_Face_Geometric_Information&& inner_face_geometeric_information)
        : cells_(cells), inner_faces_(std::move(inner_face_geometeric_information)) {};

    std::vector<Residual> calculate_RHS(const std::vector<Solution>& solutions) const {
        const auto num_solution = solutions.size();
        std::vector<Residual> RHS(num_solution, 0);
        this->inner_faces_.calculate_RHS_inner_face<Numerical_Flux_Function>(RHS, solutions);
        this->cells_.scale_RHS(RHS);
        return RHS_faces;
    }

private:
    const Cells& cells_;
    const Inner_Faces inner_faces_;
};

template<typename G, typename T, typename NF>
class Discrete_Equation
{
    using Governing_Equation = G;
    using Time_Integral_Method = T;
    using Numerical_Flux_Function = NF;
public:
    Discrete_Equation(Geometric_Information&& geometric_information)
        : Cells(std::move(geometric_information.cell_geometric_information))
        semi_discrete_equation_(this->cells_, std::move(geometric_information.inner_face_geometric_information)) {};

    template<typename End_Condition, typename Post_Condition>
    void solve(std::vector<G::Solution>& solutions) {
        // post initial solution
        while (true) {
            const auto projected_maximum_lambdas = Governing_Equation::calculate_coordinate_projected_maximum_lambdas(solutions);
            auto time_step = this->cells_.calculate_time_step(projected_maximum_lambdas);

            if (End_Condition::check(time_step)) {
                Time_Integral_Method::update_solutions(this->semi_discrete_equation_, solutions, time_step);
                //post last solution
                break;
            }

            if (Post_Condition::check(time_step)) {
                Time_Integral_Method::update_solutions(this->semi_discrete_equation_, solutions, time_step);
                //post
            }
            else
                Time_Integral_Method::update_solutions(this->semi_discrete_equation_, solutions, time_step);
        }
    }

private:
    const Cells cells_;
    const Semi_Discrete_Equation<Numerical_Flux_Function> semi_discrete_equation_;
};