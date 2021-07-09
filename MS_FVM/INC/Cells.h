#pragma once
#include "Grid_Information_Builder.h"
#include "Governing_Equation.h"

template <typename G>
class Cells
{    
    static_require(ms::is_Gov_Eq<G>, "Wrong governing equation");
    
    using Gov_Eq = G;
    using Solution = G::Solution;
    using Residual = EuclideanVector<Gov_Eq::num_equation_>;

public:
    Cells(Cell_Grid_Information&& cell_grid_information);

    double calculate_time_step(const std::vector<Solution>& solutions) const;
    void scale_RHS(std::vector<Residual>& RHS) const;

private:
    const size_t num_cell_ = 0;
    const std::vector<double> volumes_;
    const std::vector<std::array<double, Gov_Eq::physical_domain_dimension_>> coordinate_projected_volumes_;
    std::vector<double> residual_scale_factors_;
};


//template definition
template <typename G>
Cells<G>::Cells(Cell_Grid_Information&& cell_grid_information)
    : volumes_(std::move(cell_grid_information.volumes)),
    coordinate_projected_volumes_(std::move(cell_grid_information.coordinate_projected_volumes)),
    num_cell_(this->volumes_.size()) {

    this->residual_scale_factors_.resize(num_cell_);
    for (size_t i = 0; i < this->num_cell_; ++i)
        this->residual_scale_factors_[i] = -1.0 / this->volumes_[i];
};


template <typename G>
double Cells<G>::calculate_time_step(const std::vector<Solution>& solutions) const {
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
void Cells<G>::scale_RHS(std::vector<Residual>& RHS) const {
    for (size_t i = 0; i < this->num_cell_; ++i)
        RHS[i] *= this->residual_scale_factors_[i];
}