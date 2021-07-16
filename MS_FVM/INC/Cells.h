#pragma once
#include "Grid_Data_Processor.h"
#include "Governing_Equation.h"
#include "Spatial_Discrete_Method.h"
#include "Reconstruction_Method.h"

template <typename Governing_Equation, typename Spatial_Discrete_Method, typename Reconstruction_Method>
class Cells;


//Governing Equation, Reconstruction Method에 관계없이
//FVM이면 공통으로 사용하는 variable & method
template <size_t space_dimension>
class Cells_FVM_Base
{
    using SpaceVector = EuclideanVector<space_dimension>;
public:
    Cells_FVM_Base(const Processed_Grid_Data<space_dimension>& processed_grid_data);

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
    Cells_FVM_MLP_Base(const Processed_Grid_Data<space_dimension>& processed_grid_data);

protected:
    std::vector<std::vector<size_t>> vertex_node_indexes_set_;
    std::vector<std::vector<size_t>> neighbor_cell_container_indexes_set_;
    std::vector<Dynamic_Matrix_> center_to_center_matrixes_;
    std::vector<Dynamic_Matrix_> center_to_vertex_matrixes_;
    std::unordered_map<size_t, std::set<size_t>> vertex_node_index_to_neighbor_cell_container_indexes;
};



//지배방정식에 관계없이 
//FVM이고 MLPu1이면 공통으로 사용하는 method
template <typename Gradient_Method, size_t space_dimension>
class Cells_FVM_MLP_u1 : public Cells_FVM_MLP_Base<space_dimension>
{
public:
    template <typename Solution>
    std::vector<Matrix<Solution::dimension(), space_dimension>> calculate_gradient(const std::vector<Solution>& solutions) const;

private:
    template <typename Solution>
    std::vector<Dynamic_Matrix_> calculate_solution_delta_matrixes(const std::vector<Solution>& solutions) const;
    template <typename Solution>
    std::unordered_map<size_t, std::pair<Solution, Solution>> calculate_vertex_node_index_to_min_max_solution(const std::vector<Solution>& solutions) const;
    template <typename Solution>
    Dynamic_Matrix_ calculate_vertex_solution_matrix(const Solution& solution, const size_t num_vertex) const;
    double MLP_u1(const double vertex_solution_delta, const double vertex_solution, const double min_solution, const double max_solution) const;
};


template <typename Governing_Equation>
class Cells<Governing_Equation, FVM, Constant_Reconstruction> : public Cells_FVM_Base<2>
{
public:
    Cells(const Processed_Grid_Data<2>& processed_grid_data): Cells_FVM_Base<2>(processed_grid_data) {};
};


template <typename Governing_Equation, typename Gradient_Method>
class Cells<Governing_Equation, FVM, MLP_u1<Gradient_Method>> : public Cells_FVM_MLP_u1<Gradient_Method,2>
{
public:
    Cells(const Processed_Grid_Data<2>& processed_grid_data) : Cells_FVM_MLP_u1<Gradient_Method, 2>(processed_grid_data) {};
};


//template definition
template <size_t space_dimension>
Cells_FVM_Base<space_dimension>::Cells_FVM_Base(const Processed_Grid_Data<space_dimension>& processed_grid_data){
    const auto& cell_geometries = processed_grid_data.cell_geometries;
    this->num_cell_ = cell_geometries.size();

    this->centers_.reserve(this->num_cell_);
    this->volumes_.reserve(this->num_cell_);
    this->coordinate_projected_volumes_.reserve(this->num_cell_);
    this->residual_scale_factors_.reserve(this->num_cell_);
    for (const auto& geometry : cell_geometries) {
        const auto volume = geometry.volume();

        this->centers_.push_back(geometry.center_node());
        this->volumes_.push_back(volume);
        this->coordinate_projected_volumes_.push_back(geometry.coordinate_projected_volume());
        this->residual_scale_factors_.push_back(1.0 / volume);
    }
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
    std::cout << "============================================================\n";
    std::cout << "\t\t Error Anlysis\n";
    std::cout << "============================================================\n";

    if constexpr (std::is_same_v<Governing_Equation, Linear_Advection_2D> ) {
        const auto exact_solutions = Initial_Condition::template calculate_exact_solutions<Governing_Equation>(this->centers_, time);
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
    else
        std::cout << Governing_Equation::name() << " does not provide error analysis result.\n\n";
}


template <size_t space_dimension>
Cells_FVM_MLP_Base<space_dimension>::Cells_FVM_MLP_Base(const Processed_Grid_Data<space_dimension>& processed_grid_data)
    :Cells_FVM_Base<space_dimension>(processed_grid_data) {
    this->vertex_node_indexes_set_.reserve(this->num_cell_);
    this->neighbor_cell_container_indexes_set_.reserve(this->num_cell_);
    this->vertex_node_index_to_neighbor_cell_container_indexes.reserve(this->num_cell_);
    this->center_to_center_matrixes_.reserve(this->num_cell_);
    this->center_to_vertex_matrixes_.reserve(this->num_cell_);

    
    const auto& cell_geometries = processed_grid_data.cell_geometries;
    const auto& vertex_node_index_to_neighbor_cell_container_indexes = processed_grid_data.vertex_node_index_to_cell_container_indexes;
    for (size_t i = 0; i < this->num_cell_; ++i) {
        const auto& geometry = cell_geometries[i];
        const auto center_node = geometry.center_node();

        // cells_vertex_node_indexes
        auto vertex_node_indexes = geometry.vertex_node_indexes(); 

        // cells neighbor cell container indexes
        std::set<size_t> neighbor_cell_container_indexes_temp;
        for (const auto vertex_node_index : vertex_node_indexes) {
            const auto& cell_container_indexes = vertex_node_index_to_neighbor_cell_container_indexes.at(vertex_node_index);
            neighbor_cell_container_indexes_temp.insert(cell_container_indexes.begin(), cell_container_indexes.end());
        }
        neighbor_cell_container_indexes_temp.erase(i);
        std::vector<size_t> neighbor_cell_container_indexes(neighbor_cell_container_indexes_temp.begin(), neighbor_cell_container_indexes_temp.end());

        //center to center matrix
        const auto num_neighbor_cell = neighbor_cell_container_indexes.size();

        Dynamic_Matrix_ center_to_center_matrix(space_dimension, num_neighbor_cell);
        for (size_t i = 0; i < num_neighbor_cell; ++i) {
            const auto neighbor_geometry = cell_geometries[neighbor_cell_container_indexes[i]];
            const auto neighbor_center = neighbor_geometry.center_node();
            const auto center_to_center = center_node - neighbor_center;
            for (size_t j = 0; j < space_dimension; ++j)
                center_to_center_matrix.at(j, i) = center_to_center[j];
        }

        //center to vertex matrix
        const auto vertex_nodes = geometry.vertex_nodes();
        const auto num_vertex = vertex_nodes.size();

        Dynamic_Matrix_ center_to_vertex_matrix(space_dimension, num_vertex);
        for (size_t i = 0; i < num_vertex; ++i) {
            const auto center_to_vertex = center_node - vertex_nodes[i];
            for (size_t j = 0; j < space_dimension; ++j)
                center_to_vertex_matrix.at(j, i) = center_to_vertex[j];
        }
        
        this->vertex_node_indexes_set_.push_back(std::move(vertex_node_indexes));
        this->neighbor_cell_container_indexes_set_.push_back(std::move(neighbor_cell_container_indexes));
        this->center_to_center_matrixes_.push_back(std::move(center_to_center_matrix));
        this->center_to_vertex_matrixes_.push_back(std::move(center_to_vertex_matrix));
        this->vertex_node_index_to_neighbor_cell_container_indexes = vertex_node_index_to_neighbor_cell_container_indexes;
    }        
}



template <typename Gradient_Method, size_t space_dimension>
template <typename Solution>
std::vector<Matrix<Solution::dimension(), space_dimension>> Cells_FVM_MLP_u1<Gradient_Method, space_dimension>::calculate_gradient(const std::vector<Solution>& solutions) const {
    const auto solution_delta_matrixes = this->calculate_solution_delta_matrixes(solutions);    
    auto solution_gradients = Gradient_Method::solution_gradients(this->center_to_center_matrixes_, solution_delta_matrixes);

    const auto vertex_node_index_to_min_max_solution = this->calculate_vertex_node_index_to_min_max_solution(solutions);

    constexpr size_t num_equation = Solution::dimension();
    for (size_t i = 0; i < this->num_cell_; ++i) {
        const auto num_vertex = this->vertex_node_indexes_set_[i].size();
        const auto vertex_solution_matrix = this->calculate_vertex_solution_matrix(solutions[i], num_vertex);

        auto& gradient = solution_gradients[i];

        const auto vertex_solution_delta_matrix = gradient * this->center_to_vertex_matrixes_[i];

        std::array<double, num_equation> limiting_values;
        limiting_values.fill(1);

        const auto& vertex_node_indexes = this->vertex_node_indexes_set_[i];
        for (size_t i = 0; i < num_vertex; ++i) {
            const auto vertex_node_index = vertex_node_indexes[i];
            const auto& [min_solution, max_solution] = vertex_node_index_to_min_max_solution.at(vertex_node_index);

            for (size_t e = 0; e < num_equation; ++e) {
                const auto limiting_value = this->MLP_u1(vertex_solution_matrix.at(e, i), vertex_solution_delta_matrix.at(e, i), min_solution.at(e), max_solution.at(e));
                limiting_values[e] = min(limiting_values[e], limiting_value);
            }
        }
        
        for (size_t i = 0; i < num_equation; ++i)
            for (size_t j = 0; j < space_dimension; ++j)
                gradient.at(i, j) *= limiting_values.at(i);
    }


    std::vector<Matrix<num_equation, space_dimension>> limited_solution_gradients;
    limited_solution_gradients.reserve(this->num_cell_);

    for (const auto& solution_gradient : solution_gradients)
        limited_solution_gradients.push_back(solution_gradient);

    return limited_solution_gradients;
}

template <typename Gradient_Method, size_t dim>
template <typename Solution>
std::vector<Dynamic_Matrix_> Cells_FVM_MLP_u1<Gradient_Method, dim>::calculate_solution_delta_matrixes(const std::vector<Solution>& solutions) const {
    std::vector<Dynamic_Matrix_> solution_delta_matrixes;
    solution_delta_matrixes.reserve(this->num_cell_);

    const size_t num_equation = Solution::dimension();
    for (size_t i = 0; i < this->num_cell_; ++i) {
        const auto neighbor_indexes = this->neighbor_cell_container_indexes_set_.at(i);
        const auto num_neighbor = neighbor_indexes.size();

        Dynamic_Matrix_ solution_delta_matrix(num_equation, num_neighbor);
        for (size_t j = 0; j < num_neighbor; ++j) {
            const auto solution_delta = solutions[neighbor_indexes[j]] - solutions[i];
            for (size_t k = 0; k < num_equation; ++k)
                solution_delta_matrix.at(k, j) = solution_delta[k];
        }
        solution_delta_matrixes.push_back(std::move(solution_delta_matrix));
    }

    return solution_delta_matrixes;
}

template <typename Gradient_Method, size_t dim>
template <typename Solution>
std::unordered_map<size_t, std::pair<Solution, Solution>> Cells_FVM_MLP_u1<Gradient_Method, dim>::calculate_vertex_node_index_to_min_max_solution(const std::vector<Solution>& solutions) const {
    constexpr size_t num_equation = Solution::dimension();
    
    std::unordered_map<size_t, std::pair<Solution, Solution>> vertex_node_index_to_min_max_solution;
    for (const auto& [vertex_node_index, neighbor_cell_container_indexes] : this->vertex_node_index_to_neighbor_cell_container_indexes) {

        std::array<std::vector<double>, num_equation> equation_wise_solutions;
        for (const auto& cell_container_index : neighbor_cell_container_indexes) {
            for (size_t i = 0; i < num_equation; ++i)
                equation_wise_solutions[i].push_back(solutions[cell_container_index][i]);
        }

        std::array<double, num_equation> min_solution;
        std::array<double, num_equation> max_solution;
        for (size_t i = 0; i < num_equation; ++i) {
            min_solution[i] = *std::min_element(equation_wise_solutions[i].begin(), equation_wise_solutions[i].end());
            max_solution[i] = *std::max_element(equation_wise_solutions[i].begin(), equation_wise_solutions[i].end());
        }
        Solution min_sol = min_solution;
        Solution max_sol = max_solution;

        vertex_node_index_to_min_max_solution.emplace(vertex_node_index, std::make_pair(min_sol, max_sol));
    }
    return vertex_node_index_to_min_max_solution;
}


template <typename Gradient_Method, size_t dim>
template <typename Solution>
Dynamic_Matrix_ Cells_FVM_MLP_u1<Gradient_Method, dim>::calculate_vertex_solution_matrix(const Solution& solution, const size_t num_vertex) const {
    constexpr size_t num_equation = Solution::dimension();

    Dynamic_Matrix_ vertex_solution_matrix(num_equation, num_vertex);
    for (size_t i = 0; i < num_equation; ++i)
        for (size_t j = 0; j < num_vertex; ++j)
            vertex_solution_matrix.at(i, j) = solution[i];

    return vertex_solution_matrix;
}

template <typename Gradient_Method, size_t dim>
double Cells_FVM_MLP_u1<Gradient_Method, dim>::MLP_u1(const double vertex_solution_delta, const double vertex_solution, const double min_solution, const double max_solution) const {
    if (vertex_solution_delta < 0)
        return min((min_solution - vertex_solution) / vertex_solution_delta, 1);
    else
        return min((max_solution - vertex_solution) / vertex_solution_delta, 1);
}