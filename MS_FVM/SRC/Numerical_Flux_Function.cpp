#include "../INC/Numerical_Flux_Function.h"

std::vector<LLF<Euler_2D>::Numerical_Flux_> LLF<Euler_2D>::calculate(const std::vector<Solution_>& conservative_variables, const std::vector<Space_Vector_>& normals, const std::vector<std::pair<size_t, size_t>>& oc_nc_index_pairs) {
    const auto num_cell = conservative_variables.size();

    std::vector<Solution_> primitive_variables(num_cell);
    for (size_t i = 0; i < num_cell; ++i)
        primitive_variables[i] = Euler_2D::conservative_to_primitive(conservative_variables[i]);

    const auto physical_fluxes = Euler_2D::physical_fluxes(conservative_variables, primitive_variables);

    const auto num_inner_face = normals.size();

    std::vector<Numerical_Flux_> inner_face_numerical_fluxes(num_inner_face);
    for (size_t i = 0; i < num_inner_face; ++i) {
        const auto [oc_index, nc_index] = oc_nc_index_pairs[i];
        const auto& oc_physical_flux = physical_fluxes[oc_index];
        const auto& nc_physical_flux = physical_fluxes[nc_index];

        const auto& oc_side_cvariable = conservative_variables[oc_index];
        const auto& nc_side_cvariable = conservative_variables[nc_index];
        const auto& oc_side_pvariable = primitive_variables[oc_index];
        const auto& nc_side_pvariable = primitive_variables[nc_index];
        const auto& normal = normals[i];
        const auto inner_face_maximum_lambda = Euler_2D::inner_face_maximum_lambda(oc_side_pvariable, nc_side_pvariable, normal);

        inner_face_numerical_fluxes[i] = 0.5 * ((oc_physical_flux + nc_physical_flux) * normal + inner_face_maximum_lambda * (oc_side_cvariable - nc_side_cvariable));
    }
    return inner_face_numerical_fluxes;
};

LLF<Euler_2D>::Numerical_Flux_ LLF<Euler_2D>::calculate(const Solution_& oc_side_cvariable, const Solution_& nc_side_cvariable, const Space_Vector_& normal) {
    const auto oc_side_pvariable = Euler_2D::conservative_to_primitive(oc_side_cvariable);
    const auto nc_side_pvariable = Euler_2D::conservative_to_primitive(nc_side_cvariable);

    const auto oc_physical_flux = Euler_2D::physical_flux(oc_side_cvariable, oc_side_pvariable);
    const auto nc_physical_flux = Euler_2D::physical_flux(nc_side_cvariable, nc_side_pvariable);
    const auto inner_face_maximum_lambda = Euler_2D::inner_face_maximum_lambda(oc_side_pvariable, nc_side_pvariable, normal);

    Numerical_Flux_ LLF_flux = 0.5 * ((oc_physical_flux + nc_physical_flux) * normal + inner_face_maximum_lambda * (oc_side_cvariable - nc_side_cvariable));
    return LLF_flux;
}