#pragma once
#include "Grid_Info_Extractor.h"
#include "Numerical_Flux_Function.h"
#include "Spatial_Discrete_Method.h"


template <typename SDM, typename NF>
class Inner_Faces;

template<typename NF>
class Inner_Faces_FVM_Base
{
    static_require(ms::is_numeirical_flux_function<NF>, "Wrong numerical flux");

    using Numerical_Flux_Function   = NF;
    using Space_Vector              = typename Numerical_Flux_Function::Governing_Equation::Space_Vector;

    static constexpr size_t dimensions_ = Numerical_Flux_Function::Governing_Equation::dimension();

public:
    Inner_Faces_FVM_Base(Grid_Info_Extractor<FVM, dimensions_>::Inner_Face_Infos&& inner_face_grid_info);

protected:
    size_t num_inner_face_ = 0;
    const std::vector<Space_Vector> normals_;
    const std::vector<std::pair<size_t, size_t>> owner_neighbor_cell_container_indexes_;
    const std::vector<double> areas_;
};

template<typename NF>
class Inner_Faces<FVM, NF> : public Inner_Faces_FVM_Base<NF>
{
    static_require(ms::is_numeirical_flux_function<NF>, "Wrong numerical flux");

    using Numerical_Flux_Function = NF;
    using Solution = Numerical_Flux_Function::Solution;
    using Residual = Numerical_Flux_Function::Numerical_Flux;

    static constexpr size_t dimensions_ = Numerical_Flux_Function::Governing_Equation::dimension();

public:
    Inner_Faces(Grid_Info_Extractor<FVM, dimensions_>::Inner_Face_Infos&& inner_face_grid_info)
        : Inner_Faces_FVM_Base<Numerical_Flux_Function>(std::move(inner_face_grid_info)) {};

    void calculate_RHS(std::vector<Residual>& RHS, const std::vector<Solution>& solutions) const;
};


//template definition part
template <typename NF>
Inner_Faces_FVM_Base<NF>::Inner_Faces_FVM_Base(Grid_Info_Extractor<FVM, dimensions_>::Inner_Face_Infos&& inner_face_info)
    : normals_(std::move(inner_face_info.normals)), owner_neighbor_cell_container_indexes_(std::move(inner_face_info.owner_neighbor_container_indexes)), areas_(std::move(inner_face_info.areas)) {
    this->num_inner_face_ = this->normals_.size();
};

template <typename NF>
void Inner_Faces<FVM, NF>::calculate_RHS(std::vector<Residual>& RHS, const std::vector<Solution>& solutions) const {
    const auto numerical_fluxes = Numerical_Flux_Function::calculate(solutions, this->normals_, this->owner_neighbor_cell_container_indexes_);
    for (size_t i = 0; i < this->num_inner_face_; ++i) {
        const auto [container_index_o, container_index_n] = this->owner_neighbor_cell_container_indexes_[i];
        const auto delta_RHS = this->areas_[i] * numerical_fluxes[i];
        RHS[container_index_o] -= delta_RHS;
        RHS[container_index_n] += delta_RHS;
    }
};
