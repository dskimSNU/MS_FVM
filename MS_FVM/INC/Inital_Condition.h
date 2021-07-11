#pragma once
#include <numbers>
#include "EuclideanVector.h"

class Sine_Wave_2D {
    static constexpr size_t num_eqation_ = 1;
    static constexpr size_t dimension_ = 2;
    static constexpr double pi_ = std::numbers::pi;
    using Space_Vector = EuclideanVector<dimension_>;
    using Solution = EuclideanVector<num_eqation_>;
public:
    static std::vector<Solution> calculate_solutions(const std::vector<Space_Vector>& cell_centers);

private:
    Sine_Wave_2D(void) = delete;
};