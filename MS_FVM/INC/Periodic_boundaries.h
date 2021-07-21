#pragma once

#include "Periodic_Boundaries_FVM.h"
#include "Spatial_Discrete_Method.h"


template <typename Spatial_Discrete_Method, typename Reconstruction_Method, size_t space_dimension>
class Periodic_Boundaries;


template<size_t space_dimension>
class Periodic_Boundaries<FVM, Constant_Reconstruction, space_dimension> : public Periodic_Boundaries_FVM_Constant<space_dimension>
{
public:
    Periodic_Boundaries(Grid<space_dimension>&& grid) : Periodic_Boundaries_FVM_Constant<space_dimension>(std::move(grid)) {};
};


template<typename Reconstruction_Method, size_t space_dimension>
class Periodic_Boundaries<FVM, Reconstruction_Method, space_dimension> : public Periodic_Boundaries_FVM_Linear<space_dimension>
{
public:
    Periodic_Boundaries(Grid<space_dimension>&& grid) : Periodic_Boundaries_FVM_Linear<space_dimension>(std::move(grid)) {};
};