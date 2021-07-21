#pragma once
#include "Grid_Builder.h"

//FVM이면 공통으로 사용하는 variable & method
template <size_t space_dimension>
class Cells_FVM
{
    using SpaceVector = EuclideanVector<space_dimension>;
public:
    Cells_FVM(const Grid<space_dimension>& grid);

    double calculate_time_step(const std::vector<std::array<double, space_dimension>>& coordinate_projected_maximum_lambdas, const double cfl) const;

    template <typename Residual>
    void scale_RHS(std::vector<Residual>& RHS) const;

    template <typename Initial_Condtion>
    auto calculate_initial_solutions(void) const;

    template <typename Initial_Condition, typename Governing_Equation, typename Solution>
    void estimate_error(const std::vector<Solution>& computed_solution, const double time) const;


protected:
    size_t num_cell_ = 0;
    std::vector<SpaceVector> centers_;
    std::vector<double> volumes_;
    std::vector<std::array<double, space_dimension>> coordinate_projected_volumes_;
    std::vector<double> residual_scale_factors_;
};


//template definition
template <size_t space_dimension>
Cells_FVM<space_dimension>::Cells_FVM(const Grid<space_dimension>& grid) {
    SET_TIME_POINT;

    const auto& cell_elements = grid.elements.cell_elements;
    this->num_cell_ = cell_elements.size();

    this->centers_.reserve(this->num_cell_);
    this->volumes_.reserve(this->num_cell_);
    this->coordinate_projected_volumes_.reserve(this->num_cell_);
    this->residual_scale_factors_.reserve(this->num_cell_);
    for (const auto& cell_elemnt : cell_elements) {
        const auto& geometry = cell_elemnt.geometry_;

        const auto volume = geometry.volume();

        this->centers_.push_back(geometry.center_node());
        this->volumes_.push_back(volume);
        this->coordinate_projected_volumes_.push_back(geometry.coordinate_projected_volume());
        this->residual_scale_factors_.push_back(1.0 / volume);
    }

    Log::content_ << std::left << std::setw(50) << "@ Cells FVM precalculation" << " ----------- " << GET_TIME_DURATION << "s\n\n";
    Log::print();
};

template <size_t space_dimension>
double Cells_FVM<space_dimension>::calculate_time_step(const std::vector<std::array<double, space_dimension>>& coordinate_projected_maximum_lambdas, const double cfl) const {
    std::vector<double> local_time_step(this->num_cell_);
    for (size_t i = 0; i < this->num_cell_; ++i) {
        const auto [x_projected_volume, y_projected_volume] = this->coordinate_projected_volumes_[i];
        const auto [x_projeced_maximum_lambda, y_projeced_maximum_lambda] = coordinate_projected_maximum_lambdas[i];

        const auto x_radii = x_projected_volume * x_projeced_maximum_lambda;
        const auto y_radii = y_projected_volume * y_projeced_maximum_lambda;

        local_time_step[i] = cfl * this->volumes_[i] / (x_radii + y_radii);
    }

    return *std::min_element(local_time_step.begin(), local_time_step.end());
}

template <size_t dim>
template <typename Residual>
void Cells_FVM<dim>::scale_RHS(std::vector<Residual>& RHS) const {
    for (size_t i = 0; i < this->num_cell_; ++i)
        RHS[i] *= this->residual_scale_factors_[i];
}

template <size_t dim>
template <typename Initial_Condtion>
auto Cells_FVM<dim>::calculate_initial_solutions(void) const {
    return Initial_Condtion::calculate_solutions(this->centers_);
}

template <size_t dim>
template <typename Initial_Condition, typename Governing_Equation, typename Solution>
void Cells_FVM<dim>::estimate_error(const std::vector<Solution>& computed_solutions, const double time) const {
    Log::content_ << "================================================================================\n";
    Log::content_ << "\t\t\t\t Error Anlysis\n";
    Log::content_ << "================================================================================\n";

    if constexpr (std::is_same_v<Governing_Equation, Linear_Advection_2D>) {
        const auto exact_solutions = Initial_Condition::template calculate_exact_solutions<Governing_Equation>(this->centers_, time);
        double global_L1_error = 0.0;
        double global_L2_error = 0.0;
        double global_Linf_error = 0.0;
        const auto num_solutions = computed_solutions.size();
        for (size_t i = 0; i < num_solutions; ++i) {
            const auto local_error = (exact_solutions[i] - computed_solutions[i]).L1_norm();
            global_L1_error += local_error;
            global_L2_error += local_error * local_error;
            global_Linf_error = max(global_Linf_error, local_error);
        }

        global_L1_error = global_L1_error / num_solutions;
        global_L2_error = global_L2_error / num_solutions;

        global_L2_error = std::sqrt(global_L2_error);

        Log::content_ << "L1 error \t\tL2 error \t\tLinf error \n";
        Log::content_ << ms::double_to_string(global_L1_error) << "\t" << ms::double_to_string(global_L2_error) << "\t" << ms::double_to_string(global_Linf_error) << "\n\n";

    }
    else
        Log::content_ << Governing_Equation::name() << " does not provide error analysis result.\n\n";

    Log::print();
}
