#pragma once
#include "Grid_Information_Builder.h"

template<typename NF>
class Inner_Faces
{
    using Numerical_Flux_Function = NF;
    using Solution = NF::Solution;
    using Residual = NF::Numerical_Flux;

public:
    Inner_Faces(Inner_Face_Grid_Information&& inner_face_grid_info)
        : normals_(std::move(inner_face_grid_info.normals)),
        owner_neighbor_cell_container_indexes(std::move(inner_face_grid_info.owner_neighbor_container_indexes)),
        areas_(std::move(inner_face_grid_info.areas)),
        num_inner_face_(this->normals_.size()) {};

    void calculate_RHS_inner_face(std::vector<Residual>& RHS, const std::vector<Solution>& solutions) const;

private:
    const size_t num_inner_face_ = 0;
    const std::vector<Physical_Domain_Vector> normals_;
    const std::vector<std::pair<size_t, size_t>> owner_neighbor_cell_container_indexes;
    const std::vector<double> areas_;
};


template<typename NF>
void Inner_Faces<NF>::calculate_RHS_inner_face(std::vector<Residual>& RHS, const std::vector<Solution>& solutions) const {
    const auto numerical_fluxes = Numerical_Flux_Function::calculate(solutions, this->normals_, this->owner_neighbor_cell_container_indexes);
    for (size_t i = 0; i < this->num_inner_face_; ++i) {
        const auto [owner_cell_container_index, neighbor_cell_container_index] = this->owner_neighbor_cell_container_indexes[i];
        const auto delta_RHS = this->areas_[i] * numerical_fluxes[i];
        RHS[owner_cell_container_index] -= delta_RHS;
        RHS[neighbor_cell_container_index] += delta_RHS;
    }
}