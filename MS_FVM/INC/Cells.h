#pragma once
#include "Grid_Builder.h"
#include "Governing_Equation.h"
#include "Spatial_Discrete_Method.h"
#include "Reconstruction_Method.h"

//Governing Equation, Reconstruction Method에 관계없이
//FVM이면 공통으로 사용하는 variable & method
template <size_t space_dimension>
class Cells_FVM_Base
{
    using SpaceVector = EuclideanVector<space_dimension>;
public:
    Cells_FVM_Base(const Grid<space_dimension>& grid);

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


//지배방정식에 관계없이 
//FVM이고 MLP계열의 Reconstruction이면 공통으로 사용하는 variable
template <size_t space_dimension>
class Cells_FVM_MLP_Base : public Cells_FVM_Base<space_dimension>
{
public:
    Cells_FVM_MLP_Base(Grid<space_dimension>&& grid);

protected:
    //Least square Gradient variable
    std::vector<std::vector<size_t>> near_cell_indexes_set_;
    std::vector<Dynamic_Matrix_> least_square_matrixes_;

    //MLP variable
    std::vector<std::vector<size_t>> vnode_indexes_set_;
    std::vector<Dynamic_Matrix_> center_to_vertex_matrixes_;
    std::unordered_map<size_t, std::set<size_t>> vnode_index_to_share_cell_indexes_;
};


template <typename Governing_Equation, typename Spatial_Discrete_Method, typename Reconstruction_Method>
class Cells;


template <typename Governing_Equation>
class Cells<Governing_Equation, FVM, Constant_Reconstruction> : public Cells_FVM_Base<Governing_Equation::space_dimension()>
{
private:
    static constexpr size_t space_dimension_ = Governing_Equation::space_dimension();

public:
    Cells(const Grid<space_dimension_>& grid) : Cells_FVM_Base<space_dimension_>(grid) {};
};


template <typename Governing_Equation, typename Gradient_Method>
class Cells<Governing_Equation, FVM, MLP_u1<Gradient_Method>> : public Cells_FVM_MLP_Base<Governing_Equation::space_dimension()>
{
private:    
    static constexpr size_t space_dimension_    = Governing_Equation::space_dimension();
    static constexpr size_t num_equation_       = Governing_Equation::num_equation();

    using Solution_             = typename Governing_Equation::Solution_;
    using Solution_Gradient_    = Matrix<num_equation_, space_dimension_>;

public:
    Cells(Grid<space_dimension_>&& grid) : Cells_FVM_MLP_Base<space_dimension_>(std::move(grid)) {};

    std::vector<Solution_Gradient_> calculate_gradient(const std::vector<Solution_>& solutions) const;

private:
    std::vector<Dynamic_Matrix_> calculate_solution_delta_matrixes(const std::vector<Solution_>& solutions) const;
    std::unordered_map<size_t, std::pair<Solution_, Solution_>> calculate_vertex_node_index_to_min_max_solution(const std::vector<Solution_>& solutions) const;
    Dynamic_Matrix_ calculate_center_solution_matrix(const Solution_& solution, const size_t num_vertex) const;
    double MLP_u1(const double vertex_solution_delta, const double center_solution, const double min_solution, const double max_solution) const;
};


//template definition
template <size_t space_dimension>
Cells_FVM_Base<space_dimension>::Cells_FVM_Base(const Grid<space_dimension>& grid){
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

    Log::content_ << std::left << std::setw(50) << "@ Construct Cells FVM Base" << " ----------- " << GET_TIME_DURATION << "s\n\n";
    Log::print();
};

template <size_t space_dimension>
double Cells_FVM_Base<space_dimension>::calculate_time_step(const std::vector<std::array<double, space_dimension>>& coordinate_projected_maximum_lambdas, const double cfl) const{
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
void Cells_FVM_Base<dim>::scale_RHS(std::vector<Residual>& RHS) const {
    for (size_t i = 0; i < this->num_cell_; ++i)
        RHS[i] *= this->residual_scale_factors_[i];
}

template <size_t dim>
template <typename Initial_Condtion>
auto Cells_FVM_Base<dim>::calculate_initial_solutions(void) const {
    return Initial_Condtion::calculate_solutions(this->centers_);
}

template <size_t dim>
template <typename Initial_Condition, typename Governing_Equation, typename Solution>
void Cells_FVM_Base<dim>::estimate_error(const std::vector<Solution>& computed_solutions, const double time) const {
    Log::content_ << "================================================================================\n";
    Log::content_ << "\t\t\t\t Error Anlysis\n";
    Log::content_ << "================================================================================\n";

    if constexpr (std::is_same_v<Governing_Equation, Linear_Advection_2D> && std::is_same_v<Initial_Condition, Sine_Wave_2D> ) {
 /*       const auto exact_solutions = Initial_Condition::template calculate_exact_solutions<Governing_Equation>(this->centers_, time);
        double sum_error = 0.0;
        double sum_volume = 0.0;
        const auto num_solutions = computed_solutions.size();
        for (size_t i = 0; i < num_solutions; ++i) {
            const auto solution_diff = (exact_solutions[i] - computed_solutions[i]).L1_norm();
            sum_error += solution_diff * solution_diff * this->volumes_[i];
            sum_volume += this->volumes_[i];
        }

        Log::content_ << "L2 error \n";
        Log::content_ << ms::double_to_string(std::sqrt(sum_error / sum_volume)) << "\n\n";*/


        //// new ms error v1
        //const auto exact_solutions = Initial_Condition::template calculate_exact_solutions<Governing_Equation>(this->centers_, time);
        //double global_L1_error = 0.0;
        //double global_L2_error = 0.0;
        //double global_Linf_error = 0.0;
        //const auto num_solutions = computed_solutions.size();
        //for (size_t i = 0; i < num_solutions; ++i) {
        //    const auto solution_diff = (exact_solutions[i] - computed_solutions[i]).L1_norm();
        //    global_L1_error += solution_diff;
        //    global_L2_error += solution_diff * solution_diff;
        //    global_Linf_error = max(global_Linf_error, solution_diff);
        //}
        //global_L2_error = std::sqrt(global_L2_error);

        ////normalize w.r.t. num cell
        //global_L1_error = global_L1_error / num_solutions;
        //global_L2_error = global_L2_error / num_solutions;

        //Log::content_ << "L1 error \t\tL2 error \t\tLinf error \n";
        //Log::content_ << ms::double_to_string(global_L1_error) << "\t" << ms::double_to_string(global_L2_error) <<"\t" << ms::double_to_string(global_Linf_error) << "\n\n";

        //// new ms error v2
        //const auto exact_solutions = Initial_Condition::template calculate_exact_solutions<Governing_Equation>(this->centers_, time);
        //double global_L1_error = 0.0;
        //double global_L2_error = 0.0;
        //double global_Linf_error = 0.0;
        //const auto num_solutions = computed_solutions.size();
        //for (size_t i = 0; i < num_solutions; ++i) {
        //    const auto solution_diff = (exact_solutions[i] - computed_solutions[i]).L1_norm();
        //    global_L1_error += solution_diff;
        //    global_L2_error += solution_diff * solution_diff;
        //    global_Linf_error = max(global_Linf_error, solution_diff);
        //}
        //global_L2_error = std::sqrt(global_L2_error);
        //
        //Log::content_ << "L1 error \t\tL2 error \t\tLinf error \n";
        //Log::content_ << ms::double_to_string(global_L1_error) << "\t" << ms::double_to_string(global_L2_error) << "\t" << ms::double_to_string(global_Linf_error) << "\n\n";
    
        // new ms error v3
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

        //normalize w.r.t. num cell
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


template <size_t space_dimension>
Cells_FVM_MLP_Base<space_dimension>::Cells_FVM_MLP_Base(Grid<space_dimension>&& grid) : Cells_FVM_Base<space_dimension>(grid) {
    SET_TIME_POINT;

    this->vnode_index_to_share_cell_indexes_ = std::move(grid.connectivity.vnode_index_to_share_cell_indexes);

    this->vnode_indexes_set_.reserve(this->num_cell_);
    this->near_cell_indexes_set_.reserve(this->num_cell_);
    this->least_square_matrixes_.reserve(this->num_cell_);
    this->center_to_vertex_matrixes_.reserve(this->num_cell_);
    
    const auto& cell_elements = grid.elements.cell_elements;
    for (size_t i = 0; i < this->num_cell_; ++i) {
        const auto& element = cell_elements[i];
        const auto& geometry = cell_elements[i].geometry_;

        // cells_vertex_node_indexes
        auto vnode_indexes = element.vertex_node_indexes(); 

        // cells neighbor cell container indexes
        std::set<size_t> near_cell_indexes_temp;
        for (const auto vnode_index : vnode_indexes) {
            const auto& share_cell_indexes = this->vnode_index_to_share_cell_indexes_.at(vnode_index);
            near_cell_indexes_temp.insert(share_cell_indexes.begin(), share_cell_indexes.end());
        }
        near_cell_indexes_temp.erase(i);
        std::vector<size_t> near_cell_indexes(near_cell_indexes_temp.begin(), near_cell_indexes_temp.end());

        //least square matrix
        const auto num_neighbor_cell = near_cell_indexes.size();

        const auto this_center = geometry.center_node();

        Dynamic_Matrix_ center_to_center_matrix(space_dimension, num_neighbor_cell);
        for (size_t i = 0; i < num_neighbor_cell; ++i) {
            const auto& neighbor_geometry = cell_elements[near_cell_indexes[i]].geometry_;
            const auto neighbor_center = neighbor_geometry.center_node();
            const auto center_to_center = this_center - neighbor_center;
            for (size_t j = 0; j < space_dimension; ++j)
                center_to_center_matrix.at(j, i) = center_to_center[j];
        }

        auto& Rc = center_to_center_matrix;
        auto RcT = Rc.transpose();
        auto least_square_matrix = RcT * (Rc * RcT).be_inverse();

        //center to vertex matrix
        const auto vertex_nodes = geometry.vertex_nodes();
        const auto num_vertex = vertex_nodes.size();

        Dynamic_Matrix_ center_to_vertex_matrix(space_dimension, num_vertex);
        for (size_t i = 0; i < num_vertex; ++i) {
            const auto center_to_vertex = this_center - vertex_nodes[i];
            for (size_t j = 0; j < space_dimension; ++j)
                center_to_vertex_matrix.at(j, i) = center_to_vertex[j];
        }

        this->near_cell_indexes_set_.push_back(std::move(near_cell_indexes));
        this->least_square_matrixes_.push_back(std::move(least_square_matrix));

        this->vnode_indexes_set_.push_back(std::move(vnode_indexes));
        this->center_to_vertex_matrixes_.push_back(std::move(center_to_vertex_matrix));
    }        

    Log::content_ << std::left << std::setw(50) << "@ Construct Cells FVM MLP Base" << " ----------- " << GET_TIME_DURATION << "s\n\n";
    Log::print();
}



template <typename Governing_Equation, typename Gradient_Method>
std::vector<typename Cells<Governing_Equation, FVM, MLP_u1<Gradient_Method>>::Solution_Gradient_> Cells<Governing_Equation, FVM, MLP_u1<Gradient_Method>>::calculate_gradient(const std::vector<Solution_>& solutions) const {
    const auto solution_delta_matrixes = this->calculate_solution_delta_matrixes(solutions);    
    auto solution_gradients = Gradient_Method::solution_gradients(solution_delta_matrixes, this->least_square_matrixes_);

    const auto vnode_index_to_min_max_solution = this->calculate_vertex_node_index_to_min_max_solution(solutions);

    for (size_t i = 0; i < this->num_cell_; ++i) {
        auto& gradient = solution_gradients[i];
        const auto vertex_solution_delta_matrix = gradient * this->center_to_vertex_matrixes_[i];

        std::array<double, num_equation_> limiting_values;
        limiting_values.fill(1);

        const auto& vnode_indexes = this->vnode_indexes_set_[i];
        const auto num_vertex = vnode_indexes.size();

        for (size_t j = 0; j < num_vertex; ++j) {
            const auto vnode_index = vnode_indexes[j];
            const auto& [min_solution, max_solution] = vnode_index_to_min_max_solution.at(vnode_index);

            for (size_t e = 0; e < num_equation_; ++e) {
                const auto limiting_value = this->MLP_u1(vertex_solution_delta_matrix.at(e, j), solutions[i].at(e), min_solution.at(e), max_solution.at(e));
                limiting_values[e] = min(limiting_values[e], limiting_value);
            }
        }
        
        for (size_t i = 0; i < num_equation_; ++i)
            for (size_t j = 0; j < space_dimension_; ++j)
                gradient.at(i, j) *= limiting_values.at(i);
    }


    std::vector<Matrix<num_equation_, space_dimension_>> limited_solution_gradients;
    limited_solution_gradients.reserve(this->num_cell_);

    for (const auto& solution_gradient : solution_gradients)
        limited_solution_gradients.push_back(solution_gradient);

    return limited_solution_gradients;

}

template <typename Governing_Equation, typename Gradient_Method>
std::vector<Dynamic_Matrix_> Cells<Governing_Equation, FVM, MLP_u1<Gradient_Method>>::calculate_solution_delta_matrixes(const std::vector<Solution_>& solutions) const {
    std::vector<Dynamic_Matrix_> solution_delta_matrixes;
    solution_delta_matrixes.reserve(this->num_cell_);

    for (size_t i = 0; i < this->num_cell_; ++i) {
        const auto& near_cell_indexes = this->near_cell_indexes_set_.at(i);
        const auto num_near_cell = near_cell_indexes.size();

        Dynamic_Matrix_ solution_delta_matrix(num_equation_, num_near_cell);
        for (size_t j = 0; j < num_near_cell; ++j) {
            const auto solution_delta = solutions[near_cell_indexes[j]] - solutions[i];
            for (size_t k = 0; k < num_equation_; ++k)
                solution_delta_matrix.at(k, j) = solution_delta[k];
        }
        solution_delta_matrixes.push_back(std::move(solution_delta_matrix));
    }

    return solution_delta_matrixes;
}

template <typename Governing_Equation, typename Gradient_Method>
std::unordered_map<size_t, std::pair<typename Governing_Equation::Solution_, typename Governing_Equation::Solution_>> Cells<Governing_Equation, FVM, MLP_u1<Gradient_Method>>::calculate_vertex_node_index_to_min_max_solution(const std::vector<Solution_>& solutions) const {
    const size_t num_vnode = this->vnode_index_to_share_cell_indexes_.size();

    std::unordered_map<size_t, std::pair<Solution_, Solution_>> vnode_index_to_min_max_solution;
    vnode_index_to_min_max_solution.reserve(num_vnode);

    for (const auto& [vnode_index, share_cell_indexes] : this->vnode_index_to_share_cell_indexes_) {
        const size_t num_share_cell = share_cell_indexes.size();
        std::array<std::vector<double>, num_equation_> equation_wise_solutions;

        for (size_t i = 0; i < num_equation_; ++i)
            equation_wise_solutions[i].reserve(num_share_cell);

        for (const auto cell_index : share_cell_indexes) {
            for (size_t i = 0; i < num_equation_; ++i)
                equation_wise_solutions[i].push_back(solutions[cell_index][i]);
        }

        std::array<double, num_equation_> min_solution;
        std::array<double, num_equation_> max_solution;
        for (size_t i = 0; i < num_equation_; ++i) {
            min_solution[i] = *std::min_element(equation_wise_solutions[i].begin(), equation_wise_solutions[i].end());
            max_solution[i] = *std::max_element(equation_wise_solutions[i].begin(), equation_wise_solutions[i].end());
        }
        Solution_ min_sol = min_solution;
        Solution_ max_sol = max_solution;

        vnode_index_to_min_max_solution.emplace(vnode_index, std::make_pair(min_sol, max_sol));
    }

    return vnode_index_to_min_max_solution;
}


template <typename Governing_Equation, typename Gradient_Method>
Dynamic_Matrix_ Cells<Governing_Equation, FVM, MLP_u1<Gradient_Method>>::calculate_center_solution_matrix(const Solution_& solution, const size_t num_vertex) const {
    Dynamic_Matrix_ center_solution_matrix(num_equation_, num_vertex);
    for (size_t i = 0; i < num_equation_; ++i)
        for (size_t j = 0; j < num_vertex; ++j)
            center_solution_matrix.at(i, j) = solution[i];

    return center_solution_matrix;
}

template <typename Governing_Equation, typename Gradient_Method>
double Cells<Governing_Equation, FVM, MLP_u1<Gradient_Method>>::MLP_u1(const double vertex_solution_delta, const double center_solution, const double min_solution, const double max_solution) const {
    if (vertex_solution_delta < 0)
        return min((min_solution - center_solution) / vertex_solution_delta, 1);
    else
        return min((max_solution - center_solution) / vertex_solution_delta, 1);
}