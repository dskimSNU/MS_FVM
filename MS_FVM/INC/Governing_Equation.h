#pragma once
#include "Matrix.h"


class Gov_Eq {}; // Governing Equation


class SCL_2D : public Gov_Eq    // 2D Scalar Conservation Law 
{
protected:
    static constexpr size_t num_equation_ = 1;
    static constexpr size_t dimension_ = 2;

private:
    SCL_2D(void) = delete;

public:
    static constexpr size_t dimension(void) { return dimension_; };
    static constexpr size_t num_equation(void) { return num_equation_; };

    using Space_Vector  = EuclideanVector<dimension_>;
    using Solution      = EuclideanVector<num_equation_>;
    using Physical_Flux = Matrix<num_equation_, dimension_>;
};


namespace ms {
    template <typename T>
    inline constexpr bool is_governing_equation = std::is_base_of_v<Gov_Eq, T>;
    template <typename T>
    inline constexpr bool is_Scalar_Eq = std::is_base_of_v<SCL_2D, T>;
}


class Linear_Advection_2D : public SCL_2D
{
private:
    static constexpr std::array<double, dimension_> advection_speeds_ = { 1.0, 0.5 };

private:
    Linear_Advection_2D(void) = delete;

public:
    static constexpr std::array<double, dimension_> advection_speed(void) { return { 1.0, 0.5 }; };
    static Physical_Flux physical_flux(const Solution& solution);
    static std::vector<Physical_Flux> physical_fluxes(const std::vector<Solution>& solutions);
    static std::vector<std::array<double, dimension_>> coordinate_projected_maximum_lambdas(const std::vector<Solution>& solutions);
    static double inner_face_maximum_lambda(const Solution& solution_o, const Solution& solution_n, const Space_Vector& nomal_vector);
    static std::string name(void) { return "Linear_Advection_2D"; }; 
};


class Burgers_2D : public SCL_2D
{
private:
    Burgers_2D(void) = delete;

public:
    static Physical_Flux physical_flux(const Solution& solution);
    static std::vector<Physical_Flux> physical_fluxes(const std::vector<Solution>& solutions);
    static std::vector<std::array<double, dimension_>> coordinate_projected_maximum_lambdas(const std::vector<Solution>& solutions);
    static double inner_face_maximum_lambda(const Solution& solution_o, const Solution& solution_n, const Space_Vector& nomal_vector);
    static std::string name(void) { return "Burgers_2D"; };    
};