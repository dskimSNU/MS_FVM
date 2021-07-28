#pragma once
#include "Governing_Equation.h"
#include "Element.h"

template <typename Governing_Equation>
class Boundary_Flux_Function
{
protected:
	static_require(ms::is_governing_equation<Governing_Equation>, "It should be Governing Equation");

	static constexpr size_t num_equation_ = Governing_Equation::num_equation();

	using Space_Vector_		= typename Governing_Equation::Space_Vector_;
	using Solution_			= typename Governing_Equation::Solution_;
	using Boundary_Flux_	= EuclideanVector<num_equation_>;

public:
	virtual Boundary_Flux_ calculate(const Solution_& oc_cvariable, const Space_Vector_& normal) abstract;
};


template <typename Governing_Equation>
class Supersonic_Outlet_2D : public Boundary_Flux_Function<Governing_Equation>
{
protected:
	static_require(ms::is_governing_equation<Governing_Equation>, "It should be Governing Equation");

	static constexpr size_t num_equation_ = Governing_Equation::num_equation();

	using Space_Vector_ = typename Governing_Equation::Space_Vector_;
	using Solution_ = typename Governing_Equation::Solution_;
	using Boundary_Flux_ = EuclideanVector<num_equation_>;

public:
	Boundary_Flux_ calculate(const Solution_& solution, const Space_Vector_& normal) override;
};


class Slip_Wall_2D : public Boundary_Flux_Function<Euler_2D>
{
	Boundary_Flux_ calculate(const Solution_& oc_cvariable, const Space_Vector_& normal) override;
};



template <typename Governing_Equation>
class Boundary_Flux_Function_Factory
{
public:
	static std::unique_ptr<Boundary_Flux_Function<Governing_Equation>> make(const ElementType boundary_type);
};



//template definition part
template <typename Governing_Equation> 
Supersonic_Outlet_2D<Governing_Equation>::Boundary_Flux_ Supersonic_Outlet_2D<Governing_Equation>::calculate(const Solution_& solution, const Space_Vector_& normal) {
		const auto oc_physical_flux = Governing_Equation::physical_flux(solution);
		return oc_physical_flux * normal;
}


template <typename Governing_Equation>
std::unique_ptr<Boundary_Flux_Function<Governing_Equation>> Boundary_Flux_Function_Factory<Governing_Equation>::make(const ElementType boundary_type) {
	if constexpr (ms::is_SCL_2D<Governing_Equation>) {
		switch (boundary_type)
		{
		case ElementType::supersonic_outlet_2D:
			return std::make_unique<Supersonic_Outlet_2D<Governing_Equation>>();
		default:
			throw std::runtime_error("wrong element type");
			break;
		}
	}
	else {
		switch (boundary_type)
		{
		case ElementType::supersonic_outlet_2D:
			return std::make_unique<Supersonic_Outlet_2D<Governing_Equation>>();
		case ElementType::slip_wall_2D:
			return std::make_unique<Slip_Wall_2D>();
		default:
			throw std::runtime_error("wrong element type");
			break;
		}
	}
}