#pragma once
#include "Governing_Equation.h"


class Numerical_Flux_Function {};

namespace ms {
    template <typename T>
    inline constexpr bool is_Numerical_Flux_Func = std::is_base_of_v<Numerical_Flux_Function, T>;
}


template <typename G>
class LLF : public Numerical_Flux_Function
{
    static_require(ms::is_Gov_Eq<G>, "Wrong Governing Equation");

public:
    using Gov_Eq = G;
    using Physical_Domain_Vector = EuclideanVector<Gov_Eq::physical_domain_dimension_>;
    using Solution = Gov_Eq::Solution;
    using Numerical_Flux = EuclideanVector<Gov_Eq::num_equation_>;

public:
    static std::vector<Numerical_Flux> calculate(const std::vector<Solution>& solutions, const std::vector<Physical_Domain_Vector>& normals, const std::vector<std::pair<size_t, size_t>>& owner_neighbor_cell_container_indexes);
};

template <typename G>
std::vector<typename LLF<G>::Numerical_Flux> LLF<G>::calculate(const std::vector<Solution>& solutions, const std::vector<Physical_Domain_Vector>& normals, const std::vector<std::pair<size_t, size_t>>& owner_neighbor_cell_container_indexes) {
    const auto num_inner_face = normals.size();
    const auto physical_fluxes = Gov_Eq::physical_fluxes(solutions);

    std::vector<Numerical_Flux> inner_face_numerical_fluxes(num_inner_face);
    for (size_t i = 0; i < num_inner_face; ++i) {
        const auto [container_index_o, container_index_n] = owner_neighbor_cell_container_indexes[i];
        const auto physical_flux_o = physical_fluxes[container_index_o];
        const auto physical_flux_n = physical_fluxes[container_index_n];

        const auto solution_o = solutions[container_index_o];
        const auto solution_n = solutions[container_index_n];
        const auto normal = normals[i];
        const auto inner_face_maximum_lambda = Gov_Eq::inner_face_maximum_lambda(solution_o, solution_n, normal);

        //debug
        const auto central_flux = 0.5 * (physical_flux_o + physical_flux_n) * normal;
        const auto diffusion_term = 0.5 * inner_face_maximum_lambda * (solution_o - solution_n);
        const auto inner_face_numerical_flux = central_flux + diffusion_term;

        //const auto inner_face_numerical_flux = 0.5 * ((physical_flux_o + physical_flux_n) * normal + inner_face_maximum_lambda * (solution_o - solution_n));
        inner_face_numerical_fluxes[i] = inner_face_numerical_flux;
    }
    return inner_face_numerical_fluxes;
}