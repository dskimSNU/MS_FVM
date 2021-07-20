#pragma once
#include "Governing_Equation.h"
#include "Spatial_Discrete_Method.h"
#include "Reconstruction_Method.h"
#include "Cells_FVM.h"


template <typename Governing_Equation, typename Spatial_Discrete_Method, typename Reconstruction_Method>
class Cells;


template <typename Governing_Equation>
class Cells<Governing_Equation, FVM, Constant_Reconstruction> : public Cells_FVM<Governing_Equation::space_dimension()>
{
private:
    static constexpr size_t space_dimension_ = Governing_Equation::space_dimension();

public:
    Cells(const Grid<space_dimension_>& grid) : Cells_FVM<space_dimension_>(grid) {};
};


template <typename Governing_Equation, typename Reconstruction_Method>
class Cells<Governing_Equation, FVM, Reconstruction_Method> : public Cells_FVM<Governing_Equation::space_dimension()>
{
private:
    static constexpr size_t space_dimension_ = Governing_Equation::space_dimension();
    static constexpr size_t num_equation_ = Governing_Equation::num_equation();

    using Solution_ = typename Governing_Equation::Solution_;
    using Solution_Gradient_ = Matrix<num_equation_, space_dimension_>;

private:
    Reconstruction_Method reconstruction_method_;

public:
    Cells(Grid<space_dimension_>&& grid)
        : Cells_FVM<space_dimension_>(grid), Cells_Vertex_Least_Square<space_dimension_>(grid), Cells_MLP_u1<space_dimension_>(std::move(grid)) {};

    std::vector<Solution_Gradient_> calculate_gradient(const std::vector<Solution_>& solutions) const {
        return reconstruction_method_.caluclate_gradient<num_equation_>(solutions);
    }
};

//template <typename Governing_Equation, typename Gradient_Method>
//class Cells<Governing_Equation, FVM, Linear_Reconstruction<Gradient_Method>> : public Cells_FVM<Governing_Equation::space_dimension()>
//{
//private:
//    static constexpr size_t space_dimension_    = Governing_Equation::space_dimension();
//    static constexpr size_t num_equation_       = Governing_Equation::num_equation();
//
//    using Solution_             = typename Governing_Equation::Solution_;
//    using Solution_Gradient_    = Matrix<num_equation_, space_dimension_>;
//
//private:
//    Gradient_Method gradient_method;
//
//public:
//    Cells(Grid<space_dimension_>&& grid) : Cells_FVM<space_dimension_>(grid), gradient_method(grid) {};
//
//    std::vector<Solution_Gradient_> calculate_gradient(const std::vector<Solution_>& solutions) const {
//        auto solution_gradients_temp = gradient_method.calculate_solution_gradients(solutions);
//
//        //dynamic matrix to matrix
//        std::vector<Solution_Gradient_> solution_gradients;
//        solution_gradients.reserve(this->num_cell_);
//
//        for (const auto& solution_gradient : solution_gradients_temp)
//            solution_gradients.push_back(solution_gradient);
//
//        return solution_gradients;
//    };
//};


