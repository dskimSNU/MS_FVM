#pragma once
#include <numbers>
#include "EuclideanVector.h"

class Sine_Wave_2D {
    static constexpr size_t num_eqation_ = 1;
    static constexpr size_t physical_domain_dimension_ = 2;
    static constexpr double pi_ = std::numbers::pi;
    using Physical_Domain_Vector = EuclideanVector<physical_domain_dimension_>;
public:
    using Solution = EuclideanVector<num_eqation_>;
    static std::vector<Solution> calculate_solutions(const std::vector<Physical_Domain_Vector>& cell_centers);

private:
    Sine_Wave_2D(void) = delete;
};