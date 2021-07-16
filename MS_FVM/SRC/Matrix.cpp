#include "../INC/Matrix.h"

Dynamic_Matrix_::Matrix(const size_t num_row, const size_t num_column)
	:num_row_(num_row), num_column_(num_column) {
	this->value_.resize(num_row * num_column);
}

Dynamic_Matrix_ Dynamic_Matrix_::operator*(const Dynamic_Matrix_& other) const {
	return Dynamic_Matrix_(this->num_row_, other.num_column_, this->multiply_value(other));
}

double& Dynamic_Matrix_::at(const size_t row, const size_t column) {
	dynamic_require(this->is_in_range(row, column), "matrix indexes should not exceed given range");
	if (this->is_transposed())
		return this->value_[column * this->num_row_ + row];
	else
		return this->value_[row * this->num_column_ + column];
}

double Dynamic_Matrix_::at(const size_t row, const size_t column) const {
	dynamic_require(this->is_in_range(row, column), "matrix indexes should not exceed given range");
	if (this->is_transposed())
		return this->value_[column * this->num_row_ + row];
	else
		return this->value_[row * this->num_column_ + column];
}

Dynamic_Matrix_& Dynamic_Matrix_::be_transpose(void) {
	std::swap(this->num_row_, this->num_column_);

	if (this->is_transposed())
		this->transpose_type_ = CBLAS_TRANSPOSE::CblasNoTrans;
	else
		this->transpose_type_ = CBLAS_TRANSPOSE::CblasTrans;

	return *this;
}

Dynamic_Matrix_& Dynamic_Matrix_::be_inverse(void) {
	dynamic_require(this->is_square_matrix(), "invertable matrix should be square matrix");

	const auto ipiv = this->PLU_decomposition();

	const int matrix_layout = LAPACK_ROW_MAJOR;
	const lapack_int n = static_cast<int>(this->num_row_);
	const lapack_int lda = n;
	const lapack_int info = LAPACKE_dgetri(matrix_layout, n, this->value_.data(), lda, ipiv.data());

	dynamic_require(info == 0, "info should be 0 when success matrix inverse");
	//info > 0 "U is singular matrix in L-U decomposition"
	//info < 0 "fail to inverse the matrix"

	return *this;
}

Dynamic_Matrix_ Dynamic_Matrix_::transpose(void) const {
	auto result = *this;
	return result.be_transpose();
}

Dynamic_Matrix_ Dynamic_Matrix_::inverse(void) const {
	auto result = *this;
	return result.be_inverse();
}

std::pair<size_t, size_t> Dynamic_Matrix_::size(void) const {
	return std::make_pair(this->num_row_, this->num_column_);
}

bool Dynamic_Matrix_::is_square_matrix(void) const {
	return this->num_row_ == this->num_column_;
}

bool Dynamic_Matrix_::is_transposed(void) const {
	return this->transpose_type_ != CBLAS_TRANSPOSE::CblasNoTrans;
}

bool Dynamic_Matrix_::is_in_range(const size_t row, const size_t column) const {
	return row < this->num_row_&& column < this->num_column_;
}

std::vector<double> Dynamic_Matrix_::multiply_value(const Dynamic_Matrix_& other) const {
	if (this->num_column_ != other.num_row_)
		throw std::length_error("length is not matched");

	const CBLAS_LAYOUT layout = CBLAS_LAYOUT::CblasRowMajor;
	const CBLAS_TRANSPOSE transA = this->transpose_type_;
	const CBLAS_TRANSPOSE transB = other.transpose_type_;
	const MKL_INT m = static_cast<MKL_INT>(this->num_row_);
	const MKL_INT n = static_cast<MKL_INT>(other.num_column_);
	const MKL_INT k = static_cast<MKL_INT>(this->num_column_);
	const double alpha = 1;
	const MKL_INT lda = static_cast<MKL_INT>(this->leading_dimension());
	const MKL_INT ldb = static_cast<MKL_INT>(other.leading_dimension());
	const double beta = 0;
	const MKL_INT ldc = n;

	std::vector<double> new_value(this->num_row_ * other.num_column_);
	cblas_dgemm(layout, transA, transB, m, n, k, alpha, this->value_.data(), lda, other.value_.data(), ldb, beta, new_value.data(), ldc);

	return new_value;
}

size_t Dynamic_Matrix_::leading_dimension(void) const {
	// num column before OP()
	if (this->is_transposed())
		return this->num_row_;
	else
		return this->num_column_;
}

std::vector<int> Dynamic_Matrix_::PLU_decomposition(void) {
	dynamic_require(!this->is_transposed(), "matrix should not be transposed");

	const int matrix_layout = LAPACK_ROW_MAJOR;
	const lapack_int m = static_cast<int>(this->num_row_);
	const lapack_int n = static_cast<int>(this->num_column_);
	const lapack_int lda = n;
	std::vector<int> ipiv(std::min(m, n));
	lapack_int info = LAPACKE_dgetrf(matrix_layout, m, n, this->value_.data(), lda, ipiv.data());

	dynamic_require(0 <= info, "info should be greater than 0 when sucess PLU decomposition");
	return ipiv;
}