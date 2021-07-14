#pragma once
#include "Grid_Info_Extractor.h"
#include "Governing_Equation.h"

template <typename Governing_Equation, typename Spatial_Discrete_Method, typename Reconstruction_Method>
class Cells;


template <typename Governing_Equation>
class Cells_FVM_Constant_Base 
{
    static_require(ms::is_governing_equation<Governing_Equation>,       "Wrong governing equation");

protected:
    using Solution_              = typename Governing_Equation::Solution;
    using Residual_              = typename EuclideanVector<Governing_Equation::num_equation()>;
    
    static constexpr size_t dimension_ = Governing_Equation::dimension();

public:
    Cells_FVM_Constant_Base(Grid_Info_Extractor<FVM, Constant_Reconstruction, dimension_>::Cell_Infos&& cell_grid_information);
    double calculate_time_step(const std::vector<Solution_>& solutions, const double cfl) const;
    void scale_RHS(std::vector<Residual_>& RHS) const;

private:
    size_t num_cell_ = 0;
    const std::vector<double> volumes_;
    const std::vector<std::array<double, dimension_>> coordinate_projected_volumes_;
    std::vector<double> residual_scale_factors_;
};


template<typename Reconstruction_Method>
class Cells_LA_FVM_Base //LA : linear advection
{
    using Space_Vector_ = Linear_Advection_2D::Space_Vector;
    using Solution_     = Linear_Advection_2D::Solution;

    static constexpr size_t dimension_ = Linear_Advection_2D::dimension();

public:
    Cells_LA_FVM_Base(Grid_Info_Extractor<FVM, Constant_Reconstruction, dimension_>::Cell_Infos&& cell_infos)
        : cell_centers_(std::move(cell_infos.cell_centers)) {};

    template <typename Initial_Condition>
    void estimate_error(const std::vector<Solution_>& computed_solution, const double time) const;

private:
    std::vector<Space_Vector_> cell_centers_;
};


template <typename Reconstruction_Method>
class Cells<Linear_Advection_2D, FVM, Reconstruction_Method> : public Cells_FVM_Constant_Base<Linear_Advection_2D>, public Cells_LA_FVM_Base<Reconstruction_Method>
{
    using Solution_ = Linear_Advection_2D::Solution;

    static constexpr size_t dimension_ = Linear_Advection_2D::dimension();

public:
    Cells(Grid_Info_Extractor<FVM, Constant_Reconstruction, dimension_>::Cell_Infos&& cell_grid_information)
        : Cells_FVM_Constant_Base<Linear_Advection_2D>(std::move(cell_grid_information)),
          Cells_LA_FVM_Base<Reconstruction_Method>(std::move(cell_grid_information)){};
};


template <typename Governing_Equation>
class Cells<Governing_Equation, FVM, Constant_Reconstruction> : public Cells_FVM_Constant_Base<Governing_Equation>
{
    using Solution_ = Governing_Equation::Solution;

    static constexpr size_t dimension_ = Governing_Equation::dimension();

public:
    Cells(Grid_Info_Extractor<FVM, Constant_Reconstruction, dimension_>::Cell_Infos&& cell_grid_information)
        : Cells_FVM_Constant_Base<Governing_Equation>(std::move(cell_grid_information)) {};

    template <typename Initial_Condition>
    void estimate_error(const std::vector<Solution_>& computed_solution, const double time) const;
};





//template definition
template <typename Governing_Equation>
Cells_FVM_Constant_Base<Governing_Equation>::Cells_FVM_Constant_Base(Grid_Info_Extractor<FVM, Constant_Reconstruction, dimension_>::Cell_Infos&& cell_info)
    : volumes_(std::move(cell_info.volumes)),
    coordinate_projected_volumes_(std::move(cell_info.coordinate_projected_volumes)) {
    this->num_cell_ = this->volumes_.size();
    this->residual_scale_factors_.resize(this->num_cell_);
    for (size_t i = 0; i < this->num_cell_; ++i)
        this->residual_scale_factors_[i] = 1.0 / this->volumes_[i];
};

template <typename Governing_Equation>
double Cells_FVM_Constant_Base<Governing_Equation>::calculate_time_step(const std::vector<Solution_>& solutions, const double cfl) const {
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

template <typename Governing_Equation>
void Cells_FVM_Constant_Base<Governing_Equation>::scale_RHS(std::vector<Residual_>& RHS) const {
    for (size_t i = 0; i < this->num_cell_; ++i)
        RHS[i] *= this->residual_scale_factors_[i];
}


template <typename Governing_Equation>
template <typename Initial_Condition>
void Cells<Governing_Equation, FVM, Constant_Reconstruction>::estimate_error(const std::vector<Solution_>& computed_solutions, const double time) const {
    std::cout << "============================================================\n";
    std::cout << "\t\t Error Anlysis\n";
    std::cout << "============================================================\n";
    std::cout << Governing_Equation::name() << " does not provide error analysis result.\n\n";

    //if constexpr (std::is_same_v<Governing_Equation, Linear_Advection_2D>) {
        //const auto exact_solutions = Initial_Condition::template calculate_exact_solutions<Linear_Advection_2D>(this->cell_centers_, time);
        //double global_L1_error = 0.0;
        //double global_L2_error = 0.0;
        //const auto num_solutions = computed_solutions.size();
        //for (size_t i = 0; i < num_solutions; ++i) {
        //    const auto solution_diff = (exact_solutions[i] - computed_solutions[i]).L1_norm();
        //    global_L1_error += std::sqrt(solution_diff);
        //    global_L2_error += solution_diff;
        //}

        //std::cout << "L1 error \t\tL2 error \n";
        //std::cout << ms::double_to_string(global_L1_error / num_solutions) << "\t" << ms::double_to_string(global_L2_error / num_solutions) << "\n\n";
    //}
    //else
        //std::cout << Governing_Equation::name() << " does not provide error analysis result.\n\n";
}


template<typename Reconstruction_Method>
template <typename Initial_Condition>
void Cells_LA_FVM_Base<Reconstruction_Method>::estimate_error(const std::vector<Solution_>& computed_solutions, const double time) const {

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