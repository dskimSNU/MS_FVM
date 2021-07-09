#pragma once
#include "Matrix.h"

#define static_require static_assert

class Governing_Equation {};

namespace ms {
    template <typename T>
    inline constexpr bool is_Gov_Eq = std::is_base_of_v<Governing_Equation, T>;
}

class Scalar_Conservation_Law_2D : public Governing_Equation 
{
public:
    static constexpr size_t num_equation_ = 1;
    static constexpr size_t physical_domain_dimension_ = 2;
    using Solution = EuclideanVector<num_equation_>;
    using Physical_Flux = Matrix<num_equation_, physical_domain_dimension_>;
    using Physical_Domain_Vector = EuclideanVector<physical_domain_dimension_>;
};

class Linear_Advection_2D : public Scalar_Conservation_Law_2D
{
private:
    Linear_Advection_2D(void) = delete;

public:
    static std::vector<Physical_Flux> physical_fluxes(const std::vector<Solution>& solutions);
    static std::vector<std::array<double, physical_domain_dimension_>> coordinate_projected_maximum_lambdas(const std::vector<Solution>& solutions);
    static double inner_face_maximum_lambda(const Solution& solution_o, const Solution& solution_n, const Physical_Domain_Vector& nomal_vector);

private:
    static constexpr std::array<double, 2> advection_speeds_ = { 1.0, 0.5 };
};

class Burgers_2D : public Scalar_Conservation_Law_2D
{
private:
    Burgers_2D(void) = delete;

public:
    static std::vector<Physical_Flux> physical_fluxes(const std::vector<Solution>& solutions);
    static std::vector<std::array<double, physical_domain_dimension_>> coordinate_projected_maximum_lambdas(const std::vector<Solution>& solutions);
    static double inner_face_maximum_lambda(const Solution& solution_o, const Solution& solution_n, const Physical_Domain_Vector& nomal_vector);
};

