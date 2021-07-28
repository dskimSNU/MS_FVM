#pragma once
#include <numbers>
#include "EuclideanVector.h"
#include "Governing_Equation.h"

class IC {}; //Initial Condition

class Sine_Wave_2D : public IC {
    static constexpr size_t num_eqation_ = 1;
    static constexpr size_t dimension_ = 2;
    static constexpr double pi_ = std::numbers::pi;

    using Space_Vector_  = EuclideanVector<dimension_>;
    using Solution      = EuclideanVector<num_eqation_>;
public:
    static std::vector<Solution> calculate_solutions(const std::vector<Space_Vector_>& cell_centers);
    static std::string name(void) { return "Sine_Wave_2D"; };
    template <typename Governing_Equation>
    static std::vector<Solution> calculate_exact_solutions(const std::vector<Space_Vector_>& cell_centers, const double end_time);

private:
    Sine_Wave_2D(void) = delete;
};


template <>
std::vector<Sine_Wave_2D::Solution> Sine_Wave_2D::calculate_exact_solutions<Linear_Advection_2D>(const std::vector<Space_Vector_>& cell_centers, const double end_time);


class Square_Wave_2D : public IC {
    static constexpr size_t num_eqation_ = 1;
    static constexpr size_t dimension_ = 2;

    using Space_Vector_ = EuclideanVector<dimension_>;
    using Solution = EuclideanVector<num_eqation_>;
public:
    static std::vector<Solution> calculate_solutions(const std::vector<Space_Vector_>& cell_centers);
    static std::string name(void) { return "Square_Wave_2D"; };
    
    template <typename Governing_Equation>
    static std::vector<Solution> calculate_exact_solutions(const std::vector<Space_Vector_>& cell_centers, const double end_time);

private:
    Square_Wave_2D(void) = delete;
};


template <>
std::vector<Square_Wave_2D::Solution> Square_Wave_2D::calculate_exact_solutions<Linear_Advection_2D>(const std::vector<Space_Vector_>& cell_centers, const double end_time);


class Modifid_SOD_2D : public IC {
    static constexpr size_t num_eqation_ = 4;
    static constexpr size_t dimension_ = 2;

    using Space_Vector_ = EuclideanVector<dimension_>;
    using Solution_     = EuclideanVector<num_eqation_>;
public:
    static std::vector<Solution_> calculate_solutions(const std::vector<Space_Vector_>& cell_centers);
    static std::string name(void) { return "Modifid_SOD"; };

private:
    Modifid_SOD_2D(void) = delete;
};


namespace ms {
    template <typename T>
    inline constexpr bool is_initial_condition = std::is_base_of_v<IC, T>;
}