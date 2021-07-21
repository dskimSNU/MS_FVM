#pragma once
#include "Inner_Faces_FVM.h"
#include "Spatial_Discrete_Method.h"


template <typename Spatial_Discrete_Method, typename Reconstruction_Method, size_t space_dimension>
class Inner_Faces;


template <size_t space_dimension>
class Inner_Faces<FVM, Constant_Reconstruction, space_dimension> : public Inner_Faces_FVM_Constant<space_dimension>
{
public:
    Inner_Faces(Grid<space_dimension>&& grid) : Inner_Faces_FVM_Constant<space_dimension>(std::move(grid)) {};
};


template<typename Reconstruction_Method, size_t space_dimension>
class Inner_Faces<FVM, Reconstruction_Method, space_dimension> : public Inner_Faces_FVM_Linear<space_dimension>
{
public:
    Inner_Faces(Grid<space_dimension>&& grid) : Inner_Faces_FVM_Linear<space_dimension>(std::move(grid)) {};
};
