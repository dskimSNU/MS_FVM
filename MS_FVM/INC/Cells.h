#pragma once
#include "Cells_FVM.h"
#include "Spatial_Discrete_Method.h"


template <typename Spatial_Discrete_Method, size_t space_dimension>
class Cells;


template <size_t space_dimension>
class Cells<FVM, space_dimension> : public Cells_FVM<space_dimension>
{
public:
    Cells(const Grid<space_dimension>& grid) : Cells_FVM<space_dimension>(grid) {};
};