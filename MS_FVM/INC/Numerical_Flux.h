#pragma once
#include "Governing_Equation.h"

class Numerical_Flux_Function {};

template <typename G>
class LLF : public Numerical_Flux_Function
{
    using Governing_Equation = G;
    using Physical_Domain_Vector = EuclideanVector<PHYSICAL_DOMAIN_DIMENSION>;
    using Solution = G::Solution;
    using Numerical_Flux = EuclideanVector<G::num_eqation_>;
public:
  static std::vector<Numerical_Flux> calculate(const std::vector<Solution>& solutions, const std::vector<Physical_Domain_Vector>& normals, const std::vector<std::pair<size_t,size_t>>& owner_neighbor_cell_container_indexes) {
    const auto physical_fluxes = Governing_Equation::calculate_physical_fluxes(solutions);

    std::vector<Numerical_Flux> inner_face_numerical_fluxes(num_inner_face);
    for (size_t i = 0; i < num_inner_face; ++i) {
      const auto [owner_container_index,neighbor_container_index] = owner_neighbor_cell_container_indexes[i];
      const auto physical_flux_o = physical_fluxes[owner_container_index];
      const auto physical_flux_n = physical_fluxes[neighbor_container_index];

      const auto solution_o = solutions[owner_container_index];
      const auto solution_n = solutions[neihgbor_index];
      const auto normal = normals[i];
      const auto inner_face_maximum_lambda = Governing_Equation::calculate_inner_face_maximum_lambda(solution_o, solution_n, normal);

      const auto inner_face_numerical_flux = 0.5 * ((physical_flux_o + physical_flux_n) * normal + inner_face_maximum_lambda * (solution_o - solution_n));
      inner_face_numerical_fluxes[i] = inner_face_numerical_flux;
    }
    return inner_face_numerical_fluxes;
  }
};