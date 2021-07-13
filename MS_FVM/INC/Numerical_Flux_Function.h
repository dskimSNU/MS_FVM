#pragma once
#include "Governing_Equation.h"


class NFF {};    // Numerical Flux Function


namespace ms {
    template <typename T>
    inline constexpr bool is_numeirical_flux_function = std::is_base_of_v<NFF, T>;
}


template <typename G>
class LLF : public NFF  // Local Lax Fridrich method
{
    static_require(ms::is_governing_equation<G>, "Wrong Governing Equation");

public:
    using Governing_Equation    = G;
    using Space_Vector          = typename Governing_Equation::Space_Vector;
    using Solution              = typename Governing_Equation::Solution;
    using Numerical_Flux        = typename EuclideanVector<Governing_Equation::num_equation()>;

public:
    //For FVM 
    static std::vector<Numerical_Flux> calculate(const std::vector<Solution>& solutions, const std::vector<Space_Vector>& normals, const std::vector<std::pair<size_t, size_t>>& owner_neighbor_cell_container_indexes) {
        const auto num_inner_face = normals.size();
        const auto physical_fluxes = Governing_Equation::physical_fluxes(solutions); // constant solution이기 때문에 physical flux를 여러번 계산할 필요가 없다.

        std::vector<Numerical_Flux> inner_face_numerical_fluxes(num_inner_face);
        for (size_t i = 0; i < num_inner_face; ++i) {
            const auto [container_index_o, container_index_n] = owner_neighbor_cell_container_indexes[i];
            const auto physical_flux_o = physical_fluxes[container_index_o];
            const auto physical_flux_n = physical_fluxes[container_index_n];

            const auto solution_o = solutions[container_index_o];
            const auto solution_n = solutions[container_index_n];
            const auto normal = normals[i];
            const auto inner_face_maximum_lambda = Governing_Equation::inner_face_maximum_lambda(solution_o, solution_n, normal);

            inner_face_numerical_fluxes[i] = 0.5 * ((physical_flux_o + physical_flux_n) * normal + inner_face_maximum_lambda * (solution_o - solution_n));
        }
        return inner_face_numerical_fluxes;
    };

    static Numerical_Flux calculate(const Solution& solution_o, const Solution& solution_n, const Space_Vector& normal) {
        const auto physical_flux_o = Governing_Equation::physical_flux(solution_o);
        const auto physical_flux_n = Governing_Equation::physical_flux(solution_n);
        const auto inner_face_maximum_lambda = Governing_Equation::inner_face_maximum_lambda(solution_o, solution_n, normal);

        return 0.5 * ((physical_flux_o + physical_flux_n) * normal + inner_face_maximum_lambda * (solution_o - solution_n));
    }
};
