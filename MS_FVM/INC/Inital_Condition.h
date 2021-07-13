#pragma once
#include <numbers>
#include "EuclideanVector.h"
#include "Governing_Equation.h"

class Sine_Wave_2D {
    static constexpr size_t num_eqation_ = 1;
    static constexpr size_t dimension_ = 2;
    static constexpr double pi_ = std::numbers::pi;

    using Space_Vector  = EuclideanVector<dimension_>;
    using Solution      = EuclideanVector<num_eqation_>;
public:
    static std::vector<Solution> calculate_solutions(const std::vector<Space_Vector>& cell_centers);
    static std::string name(void) { return "Sine_Wave_2D"; };
    template <typename Governing_Equation>
    static std::vector<Solution> calculate_exact_solutions(const std::vector<Space_Vector>& cell_centers, const double end_time);

private:
    Sine_Wave_2D(void) = delete;
};


//template specialization
template <>
std::vector<Sine_Wave_2D::Solution> Sine_Wave_2D::calculate_exact_solutions<Linear_Advection_2D>(const std::vector<Space_Vector>& cell_centers, const double end_time);