#pragma once
#include "Governing_Equation.h"


class NFF {};    // Numerical Flux Function


template <typename Governing_Equation>
class LLF : public NFF  // Local Lax Fridrich method
{
private:
    static_require(ms::is_governing_equation<Governing_Equation>, "Wrong Governing Equation");

public:
    using Space_Vector_          = typename Governing_Equation::Space_Vector_;
    using Solution_              = typename Governing_Equation::Solution_;
    using Numerical_Flux_        = EuclideanVector<Governing_Equation::num_equation()>;

public:
    //For FVM 
    static std::vector<Numerical_Flux_> calculate(const std::vector<Solution_>& solutions, const std::vector<Space_Vector_>& normals, const std::vector<std::pair<size_t, size_t>>& oc_nc_index_pairs) {
        const auto num_inner_face = normals.size();
        const auto physical_fluxes = Governing_Equation::physical_fluxes(solutions); // constant solution이기 때문에 physical flux를 여러번 계산할 필요가 없다.

        std::vector<Numerical_Flux_> inner_face_numerical_fluxes(num_inner_face);
        for (size_t i = 0; i < num_inner_face; ++i) {
            const auto [oc_index, nc_index] = oc_nc_index_pairs[i];
            const auto oc_physical_flux = physical_fluxes[oc_index];
            const auto nc_physical_flux = physical_fluxes[nc_index];

            const auto oc_solution = solutions[oc_index];
            const auto nc_solution = solutions[nc_index];
            const auto normal = normals[i];
            const auto inner_face_maximum_lambda = Governing_Equation::inner_face_maximum_lambda(oc_solution, nc_solution, normal);

            inner_face_numerical_fluxes[i] = 0.5 * ((oc_physical_flux + nc_physical_flux) * normal + inner_face_maximum_lambda * (oc_solution - nc_solution));
        }
        return inner_face_numerical_fluxes;
    };

    static Numerical_Flux_ calculate(const Solution_& solution_o, const Solution_& solution_n, const Space_Vector_& normal) {
        const auto physical_flux_o = Governing_Equation::physical_flux(solution_o);
        const auto physical_flux_n = Governing_Equation::physical_flux(solution_n);
        const auto inner_face_maximum_lambda = Governing_Equation::inner_face_maximum_lambda(solution_o, solution_n, normal);

        return 0.5 * ((physical_flux_o + physical_flux_n) * normal + inner_face_maximum_lambda * (solution_o - solution_n));
    }
};


namespace ms {
    template <typename T>
    inline constexpr bool is_numeirical_flux_function = std::is_base_of_v<NFF, T>;
}