#pragma once
#include "Grid_Info_Extractor.h"
#include "Governing_Equation.h"

template <typename SDM, typename G>
class Cells;


template <typename G>
class Cells_FVM_Base 
{
    static_require(ms::is_governing_equation<G>, "Wrong governing equation");

protected:
    using Governing_Equation    = G;
    using Solution              = typename Governing_Equation::Solution;
    using Residual              = typename EuclideanVector<Governing_Equation::num_equation()>;
    
    static constexpr size_t dimension_ = Governing_Equation::dimension();

public:
    Cells_FVM_Base(Grid_Info_Extractor<FVM, dimension_>::Cell_Infos&& cell_grid_information);
    double calculate_time_step(const std::vector<Solution>& solutions, const double cfl) const;
    void scale_RHS(std::vector<Residual>& RHS) const;

private:
    size_t num_cell_ = 0;
    const std::vector<double> volumes_;
    const std::vector<std::array<double, dimension_>> coordinate_projected_volumes_;
    std::vector<double> residual_scale_factors_;
};


template <typename G>
class Cells_FVM_Limiter_Base : Cells_FVM_Base<G>
{
    static_require(ms::is_governing_equation<G>, "Wrong governing equation");

protected:
    using Governing_Equation     = G;
    using Space_Vector           = Governing_Equation::Space_Vector;

    static constexpr size_t dimension_ = Governing_Equation::dimension();

public:
    Cells_FVM_Limiter_Base(Grid_Info_Extractor<FVM, dimension_>::Cell_Info&& cell_grid_information);

private:
    std::vector<Space_Vector> cell_centers_;
    std::vector<std::vector<size_t>> MLP_stencils_;
};


template <typename G>
class Cells<FVM, G> : public Cells_FVM_Base<G>
{
    using Governing_Equation = G;
    using Solution           = Governing_Equation::Solution;

    static constexpr size_t dimension_ = Governing_Equation::dimension();

public:
    Cells(Grid_Info_Extractor<FVM, dimension_>::Cell_Infos&& cell_grid_information)
        : Cells_FVM_Base<Governing_Equation>(std::move(cell_grid_information)) {};

    template <typename Initial_Condition>
    void estimate_error(const std::vector<Solution>& computed_solution, const double time) const;
};


template <>
class Cells<FVM, Linear_Advection_2D> : public Cells_FVM_Base<Linear_Advection_2D>
{  
    using Space_Vector = Linear_Advection_2D::Space_Vector;    

public:
    Cells(Grid_Info_Extractor<FVM, dimension_>::Cell_Infos&& cell_grid_information)
        : Cells_FVM_Base(std::move(cell_grid_information)), cell_centers_(std::move(cell_grid_information.centers)) {};

    template <typename Initial_Condition>
    void estimate_error(const std::vector<Solution>& computed_solution, const double time) const;

private:
    std::vector<Space_Vector> cell_centers_;
};


//template definition
template <typename G>
Cells_FVM_Base<G>::Cells_FVM_Base(Grid_Info_Extractor<FVM, dimension_>::Cell_Infos&& cell_info)
    : volumes_(std::move(cell_info.volumes)),
    coordinate_projected_volumes_(std::move(cell_info.coordinate_projected_volumes)) {
    this->num_cell_ = this->volumes_.size();
    this->residual_scale_factors_.resize(this->num_cell_);
    for (size_t i = 0; i < this->num_cell_; ++i)
        this->residual_scale_factors_[i] = 1.0 / this->volumes_[i];
};

template <typename G>
double Cells_FVM_Base<G>::calculate_time_step(const std::vector<Solution>& solutions, const double cfl) const {
    const auto projected_maximum_lambdas = Governing_Equation::coordinate_projected_maximum_lambdas(solutions);

    std::vector<double> local_time_step(this->num_cell_);
    for (size_t i = 0; i < this->num_cell_; ++i) {
        const auto [x_projected_volume, y_projected_volume] = this->coordinate_projected_volumes_[i];
        const auto [x_projeced_maximum_lambda, y_projeced_maximum_lambda] = projected_maximum_lambdas[i];

        const auto x_radii = x_projected_volume * x_projeced_maximum_lambda;
        const auto y_radii = y_projected_volume * y_projeced_maximum_lambda;

        local_time_step[i] = cfl * this->volumes_[i] / (x_radii + y_radii);
    }

    return *std::min_element(local_time_step.begin(), local_time_step.end());
}

template <typename G>
void Cells_FVM_Base<G>::scale_RHS(std::vector<Residual>& RHS) const {
    for (size_t i = 0; i < this->num_cell_; ++i)
        RHS[i] *= this->residual_scale_factors_[i];
}


template <typename G>
template <typename Initial_Condition>
void Cells<FVM, G>::estimate_error(const std::vector<Solution>& computed_solutions, const double time) const {
    std::cout << "============================================================\n";
    std::cout << "\t\t Error Anlysis\n";
    std::cout << "============================================================\n";
    std::cout << Governing_Equation::name() << " does not provide error analysis result.\n\n";
}

//template <>
template <typename Initial_Condition>
void Cells<FVM, Linear_Advection_2D>::estimate_error(const std::vector<Solution>& computed_solutions, const double time) const {
    
    std::cout << "============================================================\n";
    std::cout << "\t\t Error Anlysis\n";
    std::cout << "============================================================\n";

    const auto exact_solutions = Initial_Condition::template calculate_exact_solutions<Linear_Advection_2D>(this->cell_centers_, time);    
    double global_L1_error = 0.0;
    double global_L2_error = 0.0;
    const auto num_solutions = computed_solutions.size();
    for (size_t i = 0; i < num_solutions; ++i) {
        const auto solution_diff = (exact_solutions[i] - computed_solutions[i]).L1_norm();
        global_L1_error += std::sqrt(solution_diff);
        global_L2_error += solution_diff;
    }

    std::cout << "L1 error \t\tL2 error \n";
    std::cout << ms::double_to_string(global_L1_error / num_solutions) << "\t" << ms::double_to_string(global_L2_error / num_solutions) << "\n\n";
}

//template <typename G>
//class Cells<FVM,G>
//{    
//    static_require(ms::is_Gov_Eq<G>, "Wrong governing equation");
//    
//    using Governing_Equation = G;
//    using Solution = G::Solution;
//    using Residual = EuclideanVector<Governing_Equation::num_equation_>;
//
//public:
//    Cells(Grid_Data_to_Info<Governing_Equation::dimension()>::Cell_Info&& cell_grid_information);
//    double calculate_time_step(const std::vector<Solution>& solutions) const;
//    void scale_RHS(std::vector<Residual>& RHS) const;
//
//private:
//    size_t num_cell_ = 0;
//    const std::vector<double> volumes_;
//    const std::vector<std::array<double, Governing_Equation::dimension()>> coordinate_projected_volumes_;
//    std::vector<double> residual_scale_factors_;
//};


////template definition
//template <typename G>
//Cells<FVM, G>::Cells(Grid_Data_to_Info<Governing_Equation::dimension()>::Cell_Info&& cell_info)
//    : volumes_(std::move(cell_info.volumes)),
//    coordinate_projected_volumes_(std::move(cell_info.coordinate_projected_volumes)) {
//    this->num_cell_ = this->volumes_.size();
//    this->residual_scale_factors_.resize(this->num_cell_);
//    for (size_t i = 0; i < this->num_cell_; ++i)
//        this->residual_scale_factors_[i] = 1.0 / this->volumes_[i];
//        //this->residual_scale_factors_[i] = -1.0 / this->volumes_[i];//debug
//};
//
//template <typename G>
//double Cells<FVM,G>::calculate_time_step(const std::vector<Solution>& solutions) const {
//    const auto projected_maximum_lambdas = Governing_Equation::coordinate_projected_maximum_lambdas(solutions);
//
//    std::vector<double> local_time_step(this->num_cell_);
//    for (size_t i = 0; i < this->num_cell_; ++i) {
//        const auto [x_projected_volume, y_projected_volume] = this->coordinate_projected_volumes_[i];
//        const auto [x_projeced_maximum_lambda, y_projeced_maximum_lambda] = projected_maximum_lambdas[i];
//
//        const auto x_radii = x_projected_volume * x_projeced_maximum_lambda;
//        const auto y_radii = y_projected_volume * y_projeced_maximum_lambda;
//
//        local_time_step[i] = CFL * this->volumes_[i] / (x_radii + y_radii);
//    }
//
//    return *std::min_element(local_time_step.begin(), local_time_step.end());
//}
//
//template <typename G>
//void Cells<FVM,G>::scale_RHS(std::vector<Residual>& RHS) const {
//    for (size_t i = 0; i < this->num_cell_; ++i)
//        RHS[i] *= this->residual_scale_factors_[i];
//}