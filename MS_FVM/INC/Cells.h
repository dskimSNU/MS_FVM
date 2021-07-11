#pragma once
#include "Grid_Data_to_Info.h"
#include "Governing_Equation.h"
#include "Spatial_Discrete_Method.h"

template <typename SDM, typename G>
class Cells;


template <typename G>
class Cells<FVM,G>
{    
    static_require(ms::is_Gov_Eq<G>, "Wrong governing equation");
    
    using Gov_Eq = G;
    using Solution = G::Solution;
    using Residual = EuclideanVector<Gov_Eq::num_equation_>;

public:
    Cells(Grid_Data_to_Info<Gov_Eq::dimension_>::Cell_Info&& cell_grid_information);
    double calculate_time_step(const std::vector<Solution>& solutions) const;
    void scale_RHS(std::vector<Residual>& RHS) const;

private:
    size_t num_cell_ = 0;
    const std::vector<double> volumes_;
    const std::vector<std::array<double, Gov_Eq::dimension_>> coordinate_projected_volumes_;
    std::vector<double> residual_scale_factors_;
};


//template definition
template <typename G>
Cells<FVM, G>::Cells(Grid_Data_to_Info<Gov_Eq::dimension_>::Cell_Info&& cell_info)
    : volumes_(std::move(cell_info.volumes)),
    coordinate_projected_volumes_(std::move(cell_info.coordinate_projected_volumes)) {
    this->num_cell_ = this->volumes_.size();
    this->residual_scale_factors_.resize(this->num_cell_);
    for (size_t i = 0; i < this->num_cell_; ++i)
        this->residual_scale_factors_[i] = 1.0 / this->volumes_[i];
        //this->residual_scale_factors_[i] = -1.0 / this->volumes_[i];//debug
};

template <typename G>
double Cells<FVM,G>::calculate_time_step(const std::vector<Solution>& solutions) const {
    const auto projected_maximum_lambdas = Gov_Eq::coordinate_projected_maximum_lambdas(solutions);

    std::vector<double> local_time_step(this->num_cell_);
    for (size_t i = 0; i < this->num_cell_; ++i) {
        const auto [x_projected_volume, y_projected_volume] = this->coordinate_projected_volumes_[i];
        const auto [x_projeced_maximum_lambda, y_projeced_maximum_lambda] = projected_maximum_lambdas[i];

        const auto x_radii = x_projected_volume * x_projeced_maximum_lambda;
        const auto y_radii = y_projected_volume * y_projeced_maximum_lambda;

        local_time_step[i] = CFL * this->volumes_[i] / (x_radii + y_radii);
    }

    return *std::min_element(local_time_step.begin(), local_time_step.end());
}

template <typename G>
void Cells<FVM,G>::scale_RHS(std::vector<Residual>& RHS) const {
    for (size_t i = 0; i < this->num_cell_; ++i)
        RHS[i] *= this->residual_scale_factors_[i];
}