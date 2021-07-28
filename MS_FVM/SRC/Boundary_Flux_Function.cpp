#include "../INC/Boundary_Flux_Function.h"

Slip_Wall_2D::Boundary_Flux_ Slip_Wall_2D::calculate(const Solution_& oc_cvariable, const Space_Vector_& normal) {
	const auto oc_pvariable = Euler_2D::conservative_to_primitive(oc_cvariable);
	const auto p = oc_pvariable[2];
	const auto nx = normal[0];
	const auto ny = normal[1] ;

	return { 0, p * nx, p * ny, 0 };
};


