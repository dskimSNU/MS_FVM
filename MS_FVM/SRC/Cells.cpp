#include "../INC/Cells.h"

Cells::Cells(Cell_Grid_Information&& cell_grid_information)
    : volumes_(std::move(cell_grid_information.volumes)),
    coordinate_projected_volumes_(std::move(cell_grid_information.coordinate_projected_volumes)),
    num_cell_(this->volumes_.size()) {

    this->residual_scale_factors_.resize(num_cell_);
    for (size_t i = 0; i < this->num_cell_; ++i)
        this->residual_scale_factors_[i] = -1.0 / this->volumes_[i];
};

double Cells::calculate_time_step(const std::vector<std::array<double, PHYSICAL_DOMAIN_DIMENSION>>& projected_maximum_lambdas) const {
    std::vector<double> local_time_step(this->num_cell_);
    for (size_t i = 0; i < this->num_cell_; ++i) {
        const auto [x_projected_volume, y_projected_volume] = this->coordinate_projected_volumes_[i];
        const auto [x_projeced_maximum_lambda, y_projeced_maximum_lambda] = projected_maximum_lambdas[i];

        const auto x_radii = x_projected_volume * x_projeced_maximum_lambda;
        const auto y_radii = y_projected_volume * y_projeced_maximum_lambda;

        local_time_step[i] = CFL * this->volumes_[i] / (x_radii + y_radii);
    }
}