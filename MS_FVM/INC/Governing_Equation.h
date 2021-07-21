#pragma once
#include "Matrix.h"


class Gov_Eq {}; // Governing Equation


class SCL_2D : public Gov_Eq    // 2D Scalar Conservation Law 
{
protected:
    static constexpr size_t num_equation_ = 1;
    static constexpr size_t space_dimension_ = 2;

private:
    SCL_2D(void) = delete;

public:
    using Space_Vector_  = EuclideanVector<space_dimension_>;
    using Solution_      = EuclideanVector<num_equation_>;
    using Physical_Flux_ = Matrix<num_equation_, space_dimension_>;

    static constexpr size_t space_dimension(void) { return space_dimension_; };
    static constexpr size_t num_equation(void) { return num_equation_; };
};


class Linear_Advection_2D : public SCL_2D
{
private:
    static constexpr std::array<double, space_dimension_> advection_speeds_ = { 1.0, 0.5 };

private:
    Linear_Advection_2D(void) = delete;

public:
    static constexpr std::array<double, space_dimension_> advection_speed(void) { return advection_speeds_; };
    static Physical_Flux_ physical_flux(const Solution_& solution);
    static std::vector<Physical_Flux_> physical_fluxes(const std::vector<Solution_>& solutions);
    static std::vector<std::array<double, space_dimension_>> coordinate_projected_maximum_lambdas(const std::vector<Solution_>& solutions);
    static double inner_face_maximum_lambda(const Solution_& solution_o, const Solution_& solution_n, const Space_Vector_& nomal_vector);
    static std::string name(void) { return "Linear_Advection_2D"; }; 
};


class Burgers_2D : public SCL_2D
{
private:
    Burgers_2D(void) = delete;

public:
    static Physical_Flux_ physical_flux(const Solution_& solution);
    static std::vector<Physical_Flux_> physical_fluxes(const std::vector<Solution_>& solutions);
    static std::vector<std::array<double, space_dimension_>> coordinate_projected_maximum_lambdas(const std::vector<Solution_>& solutions);
    static double inner_face_maximum_lambda(const Solution_& solution_o, const Solution_& solution_n, const Space_Vector_& nomal_vector);
    static std::string name(void) { return "Burgers_2D"; };    
};


namespace ms {
    template <typename T>
    inline constexpr bool is_governing_equation = std::is_base_of_v<Gov_Eq, T>;
    template <typename T>
    inline constexpr bool is_scalar_conservation_law = std::is_base_of_v<SCL_2D, T>;
}