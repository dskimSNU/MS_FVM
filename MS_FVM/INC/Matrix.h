#pragma once
#include "EuclideanVector.h"
#include <mkl.h>


namespace ms {
	inline constexpr size_t blas_mv_criteria = 50;
}

template<size_t num_row, size_t num_column>
class Matrix
{
public:
	Matrix(void) = default;
	Matrix(const Matrix<0, 0>& dynamic_matrix);
	template <typename... Args>
	Matrix(Args... args);

	Matrix operator+(const Matrix & A) const;
	Matrix operator*(const double scalar) const;
	EuclideanVector<num_row> operator*(const EuclideanVector<num_column>&x) const;
	bool operator==(const Matrix & A) const;

	double& at(const size_t row_index, const size_t column_index);
	double at(const size_t row_index, const size_t column_index) const;
	std::string to_string(void) const;

private:
	std::array<double, num_row * num_column> values_ = { 0 };
};


using Dynamic_Matrix_ = Matrix<0, 0>;


template<>
class Matrix<0, 0>
{
	template<size_t num_row, size_t num_column>
	friend class Matrix;
private:
	CBLAS_TRANSPOSE transpose_type_ = CBLAS_TRANSPOSE::CblasNoTrans;
	size_t num_row_ = 0;
	size_t num_column_ = 0;
	std::vector<double> value_;

public:
	Matrix(const size_t num_row, const size_t num_column);
	Matrix(const size_t num_row, const size_t num_column, std::vector<double>&& value)
		: num_row_(num_row), num_column_(num_column), value_(std::move(value)) {};

	Dynamic_Matrix_ operator*(const Dynamic_Matrix_& other) const;

	double& at(const size_t row, const size_t column);
	double at(const size_t row, const size_t column) const;
	Dynamic_Matrix_& be_transpose(void);
	Dynamic_Matrix_& be_inverse(void);
	Dynamic_Matrix_ transpose(void) const;
	Dynamic_Matrix_ inverse(void) const;
	std::pair<size_t, size_t> size(void) const;

private:
	bool is_square_matrix(void) const;
	bool is_transposed(void) const;
	bool is_in_range(const size_t irow, const size_t jcolumn) const;
	std::vector<double> multiply_value(const Dynamic_Matrix_& other) const;
	size_t leading_dimension(void) const;
	std::vector<int> PLU_decomposition(void);
};


template<size_t num_row, size_t num_column>
Matrix<num_row, num_column> operator*(const double scalar, const Matrix<num_row, num_column>& A);

template<size_t num_row, size_t num_column>
std::ostream& operator<<(std::ostream& os, const Matrix<num_row, num_column>& m);


//template definition part
template<size_t num_row, size_t num_column>
Matrix<num_row,num_column>::Matrix(const Matrix<0, 0>& dynamic_matrix) {
	dynamic_require(num_row == dynamic_matrix.num_row_ && num_column == dynamic_matrix.num_column_, "both matrix should be same size");
	std::copy(dynamic_matrix.value_.begin(), dynamic_matrix.value_.end(), this->values_.begin());
}

template<size_t num_row, size_t num_column>
template <typename... Args>
Matrix<num_row, num_column>::Matrix(Args... args) : values_{ static_cast<double>(args)... } {
	static_require(sizeof...(Args) <= num_row * num_column, "Number of arguments should be less then (num row * num column)");
	static_require(ms::are_arithmetics<Args...>, "every arguments should be arithmetics");
};

template<size_t num_row, size_t num_column>
Matrix<num_row, num_column> Matrix<num_row, num_column>::operator+(const Matrix& A) const {
	Matrix result = *this;
	for (size_t i = 0; i < num_row * num_column; ++i)
		result.values_[i] += A.values_[i];
	return result;
}

template<size_t num_row, size_t num_column>
Matrix<num_row, num_column> Matrix<num_row, num_column>::operator*(const double scalar) const {
	Matrix result = *this;
	for (size_t i = 0; i < num_row * num_column; ++i)
		result.values_[i] *= scalar;
	return result;
}

template<size_t num_row, size_t num_column>
bool Matrix<num_row, num_column>::operator==(const Matrix& A) const {
	return this->values_ == A.values_;
}

template<size_t num_row, size_t num_column>
EuclideanVector<num_row> Matrix<num_row, num_column>::operator*(const EuclideanVector<num_column>& x) const {
	std::array<double, num_row> result = { 0 };
	if constexpr (num_row * num_column < ms::blas_mv_criteria) {
		for (size_t i = 0; i < num_row; ++i)
			for (size_t j = 0; j < num_column; ++j)
				result[i] += this->values_[i * num_column + j] * x[j];
	}
	return result;
}

template<size_t num_row, size_t num_column>
double& Matrix<num_row, num_column>::at(const size_t row_index, const size_t column_index) {
	return this->values_[row_index * num_column + column_index];
}

template<size_t num_row, size_t num_column>
double Matrix<num_row, num_column>::at(const size_t row_index, const size_t column_index) const {
	return this->values_[row_index * num_column + column_index];
}

template<size_t num_row, size_t num_column>
std::string Matrix<num_row, num_column>::to_string(void) const {
	std::string result;
	for (size_t i = 0; i < num_row; ++i) {
		for (size_t j = 0; j < num_column; ++j)
			result += ms::double_to_string(this->at(i, j)) + "   \t";
		result += "\n";
	}
	return result;
}

template<size_t num_row, size_t num_column>
std::ostream& operator<<(std::ostream& os, const Matrix<num_row, num_column>& m) {
	return os << m.to_string();
}

template<size_t num_row, size_t num_column>
Matrix<num_row, num_column> operator*(const double scalar, const Matrix<num_row, num_column>& A) {
	return A * scalar;
}

//#include "MathVector.h"
//
//#include <algorithm>
//
//class RowMajorMatrix
//{
//public:
//	explicit RowMajorMatrix(void) = default;
//	explicit RowMajorMatrix(const size_t num_row, const size_t num_column);
//	explicit RowMajorMatrix(const size_t matrix_order);
//	RowMajorMatrix(const size_t num_row, const size_t num_column, const MathVector& value);
//	RowMajorMatrix(const size_t num_row, const size_t num_column, MathVector&& value);
//
//	RowMajorMatrix& operator+=(const RowMajorMatrix& other);
//	RowMajorMatrix& operator*=(const RowMajorMatrix& other);
//	RowMajorMatrix operator*(const RowMajorMatrix& other) const;
//	MathVector operator*(const MathVector& x) const;
//	bool operator==(const RowMajorMatrix& other) const;
//
	//double& at(const size_t row, const size_t column);
//	double at(const size_t row, const size_t column) const;
//	RowMajorMatrix& be_inverse(void);
//	RowMajorMatrix& be_transpose(void);
//	bool compare_with_finitie_precision(const RowMajorMatrix& other, const size_t ULP_precision = 4) const;
//	RowMajorMatrix& change_column(const size_t column_index, const MathVector& value);
//	RowMajorMatrix inverse(void) const;
//	RowMajorMatrix part(const size_t row_start_index, const size_t row_end_index, const size_t column_start_index, const size_t column_end_index) const;
//	MathVector row(const size_t row_index) const;
//	std::pair<size_t, size_t> size(void) const;	
//	std::string to_string(void) const;
//	RowMajorMatrix transpose(void) const;
//
//private:
//	size_t leading_dimension(void) const;
//	MathVector multiply_value(const RowMajorMatrix& other) const;
//	std::vector<int> PLU_decomposition(void);
//	RowMajorMatrix& transpose_value(void);
//
//	bool is_square_matrix(void) const;
//	bool is_transposed(void) const;
//	void inspect_range(const size_t irow, const size_t jcolumn) const;
//	void inspect_value_size(void) const;	
//
//private:
//	CBLAS_TRANSPOSE transpose_type_ = CBLAS_TRANSPOSE::CblasNoTrans;
//	size_t num_row_ = 0;
//	size_t num_column_ = 0;
//	MathVector value_vector_;
//};
//std::ostream& operator<<(std::ostream& os, const RowMajorMatrix& x);
//
//
//
//template <typename T>
//class JacobianFunction
//{
//private:
//	std::vector<VectorFunction<T>> gradient_set_;
//
//public:
//	JacobianFunction(void) = default;
//	JacobianFunction(const VectorFunction<T>& vector_function);
//	JacobianFunction(const VectorFunction<T>& vector_function, const size_t domain_dimension);
//	JacobianFunction(std::initializer_list<VectorFunction<T>> list) : gradient_set_(list) {};
//
//	RowMajorMatrix operator()(const MathVector& variable_vector) const;
//	bool operator==(const JacobianFunction& other) const;
//
//	std::pair<size_t, size_t> size(void) const;
//	std::string to_string(void) const;
//};
//
//
//template <typename T>
//std::ostream& operator<<(std::ostream& os, const JacobianFunction<T>& Jacobian_matrix);
//
//namespace ms {
//	template <typename T>
//	MathVector Newton_Raphson(const VectorFunction<T>& vector_function, const MathVector& initial_guess, const double convergence_criteria = 1.0E-10);
//}
//
//
//
////template definition part
//template <typename T>
//JacobianFunction<T>::JacobianFunction(const VectorFunction<T>& vector_function) {
//	*this = JacobianFunction(vector_function, vector_function.domain_dimension());
//}
//
//template <typename T>
//JacobianFunction<T>::JacobianFunction(const VectorFunction<T>& vector_function, const size_t domain_dimension) {
//	const auto range_dimension = vector_function.size();
//	this->gradient_set_.reserve(range_dimension);
//
//	for (const auto& Function : vector_function)
//		gradient_set_.push_back(Function.gradient(domain_dimension));
//}
//
//template <typename T>
//RowMajorMatrix JacobianFunction<T>::operator()(const MathVector& variable_vector) const {
//	const auto [num_row, num_column] = this->size();
//
//	RowMajorMatrix value(num_row, num_column);
//	for (size_t i = 0; i < num_row; ++i)
//		for (size_t j = 0; j < num_column; ++j)
//			value.at(i, j) = this->gradient_set_[i][j](variable_vector);
//
//	return value;
//}
//
//template <typename T>
//bool JacobianFunction<T>::operator==(const JacobianFunction& other) const {
//	if (this->size() != other.size())
//		return false;
//	else
//		return this->gradient_set_ == other.gradient_set_;
//}
//
//template <typename T>
//std::pair<size_t, size_t> JacobianFunction<T>::size(void) const {
//	return { this->gradient_set_.size(), this->gradient_set_.front().size() };
//}
//
//template <typename T>
//std::string JacobianFunction<T>::to_string(void) const {
//	std::string str;
//
//	const auto [num_row, num_column] = this->size();
//	for (size_t i = 0; i < num_row; ++i)
//		for (size_t j = 0; j < num_column; ++j)
//			str += "[" + std::to_string(i) + "," + std::to_string(j) + "]  :  " + this->gradient_set_[i][j].to_string() + "\n";
//
//	return str;
//}
//
//
//template <typename T>
//std::ostream& operator<<(std::ostream& os, const JacobianFunction<T>& Jacobian_matrix) {
//	return os << Jacobian_matrix.to_string();
//}
//
//
//namespace ms {
//	template <typename T>
//	MathVector Newton_Raphson(const VectorFunction<T>& vector_function, const MathVector& initial_guess, const double convergence_criteria) {
//		// solve vector function = 0
//		JacobianFunction<T> Jacobian(vector_function);
//
//		auto solution = initial_guess;
//		for (;;) {
//			const auto solution_delta = Jacobian(solution).inverse() * vector_function(solution);
//			solution -= solution_delta;
//
//			if (solution_delta.L2_Norm() <= convergence_criteria)
//				return solution;
//		}
//	}
//}
//
//







//
//#include "MathVector.h"
//
//#include <array>
//#include <algorithm> //min max
//#include <stdexcept>
//
//
//enum class MatrixType
//{
//	Full,
//	FullUpperTriangle,
//	FullLowerTriangle,
//	PackedUpperTriangle,
//	PackedLowerTriangle
//};
//
//
//class RowMajorMatrix
//{
//private:
//	CBLAS_TRANSPOSE transpose_type_ = CBLAS_TRANSPOSE::CblasNoTrans;
//	mutable MatrixType matrix_type_ = MatrixType::Full;
	//size_t row_ = 0;
	//size_t column_ = 0;
	//MathVector value_set_;
//
//public:
	//explicit RowMajorMatrix(void) = default;
	//explicit RowMajorMatrix(const MatrixType matrix_type, const size_t num_row, const size_t num_column);
	//explicit RowMajorMatrix(const MatrixType matrix_type, const size_t matrix_order);
	//RowMajorMatrix(const MatrixType matrix_type, const size_t matrix_order, const MathVector& value);
	//RowMajorMatrix(const MatrixType matrix_type, const size_t matrix_order, MathVector&& value);
	//RowMajorMatrix(const MatrixType matrix_type, const size_t num_row, const size_t num_column, const MathVector& value);
	//RowMajorMatrix(const MatrixType matrix_type, const size_t num_row, const size_t num_column, MathVector&& value);
//
//
//	RowMajorMatrix& operator+=(const RowMajorMatrix& other);
//	RowMajorMatrix& operator-=(const RowMajorMatrix& B);
//	RowMajorMatrix& operator*=(const double scalar);
//	RowMajorMatrix& operator*=(const RowMajorMatrix& other);
//	RowMajorMatrix operator+(const RowMajorMatrix& other) const;
//	RowMajorMatrix operator-(const RowMajorMatrix& other) const;
//	RowMajorMatrix operator*(const double scalar) const;	
//	MathVector operator*(const MathVector& vec) const;
//	RowMajorMatrix operator*(const RowMajorMatrix& other) const;
//	bool operator==(const RowMajorMatrix& B) const;
//	bool operator!=(const RowMajorMatrix& B) const {
//		return !(*this == B);
//	};
//
//	
//	RowMajorMatrix& abs(void) {
//		Math::abs(this->value_set_);
//		return *this;
//	}
//
	//double& at(const size_t row, const size_t column);

	//double at(const size_t row, const size_t column) const;
//
//	RowMajorMatrix cofactor_matrix(void) const;
//
//	void change_Column(const size_t column_index, const MathVector& vector);
//
//	void change_Row(const size_t row_index, const MathVector& vector);
//
//	MathVector column(const size_t column_index) const;
//
//	RowMajorMatrix column(const size_t start_column, const size_t end_column) const;
//
//	double determinant(void);
//
//	double determinant(void) const {
//		RowMajorMatrix A = *this;
//		return A.determinant();
//	}
//
//	RowMajorMatrix& inverse(void);
//
//	RowMajorMatrix inverse(void) const {
//		RowMajorMatrix A = *this;
//		return A.inverse();
//	}
//
//	double L2_Norm(void) const {
//		return this->value_set_.L2_Norm();
//	}
//
//	double minor(const size_t deleting_row_index, const size_t deleting_column_index) const;
//
//	void Moore_Penrose_Inverse(void);
//
//	RowMajorMatrix part(const size_t start_row, const size_t end_row, const size_t start_column, const size_t end_column) const;
//
//	RowMajorMatrix part(const size_t part_matrix_order) const;
//
//	RowMajorMatrix part(const size_t deleting_row_index, const size_t deleting_column_index) const;
//
//	std::vector<int> PLU_decomposition(void);
//
//	MathVector row(const size_t row_index) const;
//
//	RowMajorMatrix row(const size_t start_row, const size_t end_row) const;
//
//	std::array<RowMajorMatrix, 3> singular_value_decomposition(void);
//
//	std::pair<size_t, size_t> size(void) const {
//		return { this->row_, this->column_ };
//	};
//
//	RowMajorMatrix& scale_Row(const size_t start_row, const size_t end_row, const double scale_factor);
//
//	RowMajorMatrix& scale_Column(const size_t start_column, const size_t end_column, const double scale_factor);
//
//	std::string to_String(void) const;
//
//	RowMajorMatrix& transpose(void);
//
//	RowMajorMatrix transpose(void) const;
//
//private:
//	CBLAS_DIAG CBALS_Diagonal_Type(void) const;
//
//	CBLAS_UPLO CBLAS_Triangle_Type(void) const;
//
//	double* data(void) {
//		return this->value_set_.data();
//	};
//
//	const double* data(void) const {
//		return this->value_set_.data();
//	};
//
//	double determinant_Full_Matrix(void);
//
//	double determinant_Full_Matrix(void) const;
//
//	double determinant_Triangle_Matrix(void) const;
//
//	void inspect_Range(const size_t irow, const size_t jcolumn) const;
//
//	void inspect_Size(void) const;
//
//	void inspect_Value(void) const;
//
//	RowMajorMatrix& inverse_Full_Matrix(void);
//
//	RowMajorMatrix& inverse_Full_Triangle_Matrix(void);
//
//	RowMajorMatrix& inverse_Packed_Triangle_Matrix(void);
//
//	bool is_Full_Storage_Type(void) const;
//
//	bool is_Full_Storage_Value(void) const;
//
//	bool is_Lower_Triangle_Value(void) const;
//
//	bool is_Lower_Part(const size_t irow, const size_t jcolumn) const;
//		
//	bool is_Packed_Storage_Value(void) const;
//
//	bool is_Square_Matrix(void) const;
//
//	bool is_transposed(void) const;
//
//	bool is_Upper_Triangle_Value(void) const;
//
//	bool is_Upper_Triangle_Type(void) const;
//
//	bool is_Unit_Diagonal(void) const;
//	
//	bool is_Upper_Part(const size_t irow, const size_t jcolumn) const;
//
//	char LAPACK_Diagonal_Type(void) const;
//
//	char LAPACK_Triangle_Type(void) const;
//
//	size_t Leading_Dimension(void) const;
//
//	RowMajorMatrix& multiply_Full_Matrix_Full_Matrix(const RowMajorMatrix& B);
//
//	RowMajorMatrix& multiply_Full_Matrix_Full_Triangle_Matrix(const RowMajorMatrix& B);
//
//	MathVector multiply_Full_Matrix_Vector(const MathVector& vec) const;
//	
//	MathVector multiply_Full_Triangle_Matrix_Vector(const MathVector& vec) const;
//
//	MathVector multiply_Packed_Triangle_Matrix_Vector(const MathVector& vec) const;
//
//	size_t order(void) const;
//
//	RowMajorMatrix& transpose_Value(void);
//
//	RowMajorMatrix transpose_Value(void) const;
//
//	RowMajorMatrix& to_Full_Matrix(void);
//
//	RowMajorMatrix to_Full_Matrix(void) const;
//
//	const RowMajorMatrix& update_Type(void) const;
//};
//
//
//inline std::ostream& operator<<(std::ostream& os, const RowMajorMatrix& A) {
//	return os << A.to_String();
//}
//
//
//inline RowMajorMatrix operator*(const double scalar, const RowMajorMatrix& A) {
//	return A * scalar;
//}
//
//
//template <typename T>
//class JacobianMatrix
//{
//private:
//	std::vector<std::vector<T>> function_set_;
//
//public:
//	JacobianMatrix(void) = default;
//
//	JacobianMatrix(const size_t num_function, const size_t num_variable) {
//		function_set_.resize(num_function);
//		for (auto& set : function_set_)
//			set.resize(num_variable);
//	}
//
//
//	RowMajorMatrix operator()(const MathVector& variable_vector) const {
//		const auto [num_row, num_column] = this->size();
//
//		RowMajorMatrix value(MatrixType::Full, num_row, num_column);
//		for (size_t i = 0; i < num_row; ++i)
//			for (size_t j = 0; j < num_column; ++j)
//				value.at(i, j) = function_set_[i][j](variable_vector);
//
//		return value;
//	}
//
//	T& at(const size_t i_index, const size_t j_index) {
//		return function_set_[i_index][j_index];
//	}
//
//	const T& at(const size_t i_index, const size_t j_index) const {
//		return function_set_[i_index][j_index];
//	}
//
//	std::pair<size_t, size_t> size(void) const {
//		return { this->function_set_.size(), this->function_set_.front().size() };
//	}
//
//	std::string to_String(void) const {
//		std::string str;
//
//		const auto [num_row, num_column] = this->size();
//		for (size_t i = 0; i < num_row; ++i)
//			for (size_t j = 0; j < num_column; ++j)
//				str += "[" + Editor::to_String(i) + "," + Editor::to_String(j) + "]  :  " + Editor::to_String(this->at(i, j)) + "\n";
//
//		return str;
//	}
//};
//
//
//template <typename T>
//inline std::ostream& operator<<(std::ostream& os, const JacobianMatrix<T>& Jacobian_matrix) {
//	return os << Jacobian_matrix.to_String();
//}
//
//namespace Math {
//	template <typename T>
	//JacobianMatrix<T> Jacobian(const VectorFunction<T>& vector_function, const size_t num_variable) {
	//	const auto num_function = vector_function.size();

	//	JacobianMatrix<T> jacobian_matrix(num_function, num_variable);
	//	for (size_t i = 0; i < num_function; ++i)
	//		for (size_t j = 0; j < num_variable; ++j)
	//			jacobian_matrix.at(i, j) = Math::differentiate(vector_function.at(i), j);

	//	return jacobian_matrix;
//	}
//}