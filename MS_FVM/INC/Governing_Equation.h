#pragma once
#include "Matrix.h"
#include "../RSC/Setting.h"

class Governing_Equation {};

class Scalar_Conservation_Law_2D : public Governing_Equation 
{
protected:
    static constexpr size_t num_equation_ = 1;
    static constexpr size_t physical_domain_dimension_ = 2;
public:
    using Solution = EuclideanVector<num_equation_>;
    using Physical_Flux = Matrix<num_equation_, physical_domain_dimension_>;
    using Physical_Domain_Vector = EuclideanVector<physical_domain_dimension_>;
};

class Linear_Advection_2D : public Scalar_Conservation_Law_2D
{
public:
    static std::vector<Physical_Flux> calculate_physical_fluxes(const std::vector<Solution>& solutions);
    static std::vector<std::array<double, physical_domain_dimension_>> calculate_coordinate_projected_maximum_lambdas(const std::vector<Solution>& solutions);
    static double calculate_inner_face_maximum_lambda(const Solution& solution_o, const Solution& solution_n, const Physical_Domain_Vector& nomal_vector);

private:
    Linear_Advection_2D(void) = delete;

    static constexpr std::array<double, 2> advection_speeds_ = { 1.0, 0.5 };
};

class Burgers_2D : public Scalar_Conservation_Law_2D
{
public:
    static std::vector<Physical_Flux> calculate_physical_fluxes(const std::vector<Solution>& solutions);
    static std::vector<std::array<double, physical_domain_dimension_>> calculate_coordinate_projected_maximum_lambdas(const std::vector<Solution>& solutions);
    static double calculate_inner_face_maximum_lambda(const Solution& solution_o, const Solution& solution_n, const Physical_Domain_Vector& nomal_vector);

private:
    Burgers_2D(void) = delete;
};