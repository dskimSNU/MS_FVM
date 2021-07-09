#pragma once

#include "Grid_Information_Builder.h"

class Cells
{    
public:
    Cells(Cell_Grid_Information&& cell_grid_information);

    double calculate_time_step(const std::vector<std::array<double, PHYSICAL_DOMAIN_DIMENSION>>& projected_maximum_lambdas) const;
    template <typename Residual>
    void scale_RHS(std::vector<Residual>& RHS) const;

private:
    const size_t num_cell_ = 0;
    const std::vector<double> volumes_;
    const std::vector<std::array<double, PHYSICAL_DOMAIN_DIMENSION>> coordinate_projected_volumes_;
    std::vector<double> residual_scale_factors_;
};


//template definition
template <typename Residual>
void Cells::scale_RHS(std::vector<Residual>& RHS) const {
    for (size_t i = 0; i < this->num_cell_; ++i)
        RHS[i] *= this->residual_scale_factors_[i]);
}
