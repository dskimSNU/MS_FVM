#pragma once
#include "Grid_Data_to_Info.h"
#include "Numerical_Flux_Function.h"
#include "Spatial_Discrete_Method.h"


template <typename SDM, typename NF>
class Inner_Faces;


template<typename NF>
class Inner_Faces<FVM, NF>
{
    static_require(ms::is_Numerical_Flux_Func<NF>, "Wrong numerical flux");
        
    using Numerical_Flux_Func = NF;
    using Space_Vector = Numerical_Flux_Func::Gov_Eq::Space_Vector;
    using Solution = Numerical_Flux_Func::Solution;
    using Residual = Numerical_Flux_Func::Numerical_Flux;

public:
    Inner_Faces(Grid_Data_to_Info<Numerical_Flux_Func::Gov_Eq::dimension_>::Inner_Face_Info&& inner_face_grid_info);

    void calculate_RHS(std::vector<Residual>& RHS, const std::vector<Solution>& solutions) const;
     
private:
    size_t num_inner_face_ = 0;
    const std::vector<Space_Vector> normals_;
    const std::vector<std::pair<size_t, size_t>> owner_neighbor_cell_container_indexes_;
    const std::vector<double> areas_;
};


//template definition part
template <typename NF>
Inner_Faces<FVM, NF>::Inner_Faces(Grid_Data_to_Info<Numerical_Flux_Func::Gov_Eq::dimension_>::Inner_Face_Info&& inner_face_info)
    : normals_(std::move(inner_face_info.normals)), owner_neighbor_cell_container_indexes_(std::move(inner_face_info.owner_neighbor_container_indexes)), areas_(std::move(inner_face_info.areas)) {
    this->num_inner_face_ = this->normals_.size();
};

template <typename NF>
void Inner_Faces<FVM, NF>::calculate_RHS(std::vector<Residual>& RHS, const std::vector<Solution>& solutions) const {
    const auto numerical_fluxes = Numerical_Flux_Func::calculate(solutions, this->normals_, this->owner_neighbor_cell_container_indexes_);
    for (size_t i = 0; i < this->num_inner_face_; ++i) {
        const auto [container_index_o, container_index_n] = this->owner_neighbor_cell_container_indexes_[i];
        const auto delta_RHS = this->areas_[i] * numerical_fluxes[i];
        RHS[container_index_o] -= delta_RHS;
        RHS[container_index_n] += delta_RHS;
    }
};
