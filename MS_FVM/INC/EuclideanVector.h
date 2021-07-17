#pragma once
#include <array>
#include <iomanip>
#include <string>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <vector>

#define static_require static_assert
#define dynamic_require(requirement, state) if (!(requirement)) throw std::runtime_error(state)

namespace ms {
	template <typename... Args>
	inline constexpr bool are_arithmetics = (... && std::is_arithmetic_v<Args>);
}

// CLASS TEMPLATE euclidean_vector
template <size_t dim>
class EuclideanVector
{
public:
	EuclideanVector(void) = default;
	EuclideanVector(const std::array<double, dim>& other) : small_buffer_(other) {};
	EuclideanVector(const std::vector<double>& other) : elements_(other) {};
	template <typename... Args>
	EuclideanVector(Args... args);

	EuclideanVector& operator+=(const EuclideanVector& y);
	EuclideanVector& operator-=(const EuclideanVector& y);
	EuclideanVector& operator*=(const double scalar);
	EuclideanVector operator+(const EuclideanVector& y) const;	
	EuclideanVector operator-(const EuclideanVector& y) const;
	EuclideanVector operator*(const double scalar) const;
	//bool operator<(const EuclideanVector& y) const;
	bool operator==(const EuclideanVector& y) const;
	double operator[](const size_t position) const;

	double at(const size_t position) const;
	EuclideanVector& be_normalize(void);
	static constexpr size_t dimension(void);
	const double* data(void) const;
	double inner_product(const EuclideanVector& y) const;
	double L1_norm(void) const;
	double norm(void) const;
	bool is_axis_translation(const EuclideanVector& other, const size_t axis_tag) const;
	std::string to_string(void) const;

private:
	std::array<double, dim> small_buffer_ = { 0 };
	std::vector<double> elements_;
};


//user-defined deduction guides
template <typename... Args>
EuclideanVector(Args... args)->EuclideanVector<sizeof...(Args)>;
EuclideanVector(const std::vector<double>& vec)->EuclideanVector<0>;


template <size_t dim> 
std::ostream& operator<<(std::ostream& os, const EuclideanVector<dim>& x);

template <size_t dim>
EuclideanVector<dim> operator*(const double constant, const EuclideanVector<dim>& x);


namespace ms {
	inline std::string double_to_string(const double val) {
		constexpr size_t precision = 16;
		std::stringstream stream;
		stream << std::setprecision(precision) << std::noshowpoint << val;
		return stream.str();
	}

	inline bool compare_double(const double d1, const double d2, const size_t ULP_precision = 4) {
		const auto lower_ULP = d1 - std::nextafter(d1, std::numeric_limits<double>::lowest());
		const auto upper_ULP = std::nextafter(d1, std::numeric_limits<double>::max()) - d1;

		return d1 - ULP_precision * lower_ULP <= d2 && d2 <= d1 + ULP_precision * upper_ULP;
	}
}


// Template Definition Part
template <size_t dim>
template <typename... Args>
EuclideanVector<dim>::EuclideanVector(Args... args) : small_buffer_{ static_cast<double>(args)... } {
	static_require(sizeof...(Args) <= dim, "Number of arguments can not exceed dimension");
	static_require(ms::are_arithmetics<Args...>, "every arguments should be arithmetics");
};


template <size_t dim>
EuclideanVector<dim>& EuclideanVector<dim>::operator+=(const EuclideanVector& y) {
	for (size_t i = 0; i < dim; ++i)
		this->small_buffer_[i] += y.small_buffer_[i];
	return *this;
}

template <size_t dim>
EuclideanVector<dim>& EuclideanVector<dim>::operator-=(const EuclideanVector& y) {
	for (size_t i = 0; i < dim; ++i)
		this->small_buffer_[i] -= y.small_buffer_[i];
	return *this;
}

template <size_t dim>
EuclideanVector<dim>& EuclideanVector<dim>::operator*=(const double scalar) {
	for (size_t i = 0; i < dim; ++i)
		this->small_buffer_[i] *= scalar;
	return *this;
}

template <size_t dim> 
EuclideanVector<dim> EuclideanVector<dim>::operator+(const EuclideanVector& y) const {
	auto result = *this;
	return result += y;
}

template <size_t dim> EuclideanVector<dim> EuclideanVector<dim>::operator-(const EuclideanVector& y) const {
	auto result = *this;
	return result -= y;
}

template <size_t dim> EuclideanVector<dim> EuclideanVector<dim>::operator*(const double scalar) const {
	auto result = *this;
	return result *= scalar;
}

//template <size_t dim>
//bool EuclideanVector<dim>::operator<(const EuclideanVector& y) const {
//	for (size_t i = 0; i < dim; ++i) {
//		if (this->small_buffer_[i] == y.small_buffer_[i])
//			continue;
//		return this->small_buffer_[i] < y.small_buffer_[i];
//	}
//	return false;
//}

template <size_t dim> 
bool EuclideanVector<dim>::operator==(const EuclideanVector& y) const {
	for (size_t i = 0; i < dim; ++i) {
		if (this->small_buffer_[i] != y.small_buffer_[i])
			return false;
	}
	return true;
}

template <size_t dim> 
double EuclideanVector<dim>::operator[](const size_t position) const {
	dynamic_require(position <= dim, "Position should be less than dimension");
	return this->small_buffer_[position];
}

template <size_t dim>
double EuclideanVector<dim>::at(const size_t position) const {
	dynamic_require(position <= dim, "Position should be less than dimension");
	return this->small_buffer_[position];
}

template <size_t dim>
EuclideanVector<dim>& EuclideanVector<dim>::be_normalize(void) {
	const auto scale_factor = 1.0 / this->norm();
	return (*this) *= scale_factor;
}

template <size_t dim>
constexpr size_t EuclideanVector<dim>::dimension(void) {
	return dim;
}


template <size_t dim>
const double* EuclideanVector<dim>::data(void) const {
	if constexpr (dim != 0)
		return this->small_buffer_.data();
	else
		return this->elements_.data();
}

template <size_t dim>
double EuclideanVector<dim>::inner_product(const EuclideanVector& y) const {
	double result = 0;
	for (size_t i = 0; i < dim; ++i)
		result += this->small_buffer_[i] * y.small_buffer_[i];
	return result;
}

template <size_t dim>
double EuclideanVector<dim>::L1_norm(void) const {
	double L1_norm = 0.0;
	for (size_t i = 0; i < dim; ++i)
		L1_norm += std::abs(this->small_buffer_[i]);
	return L1_norm;
}

template <size_t dim>
double EuclideanVector<dim>::norm(void) const {
	return std::sqrt(this->inner_product(*this));
}

template <size_t dim>
bool EuclideanVector<dim>::is_axis_translation(const EuclideanVector& other, const size_t axis_tag) const {
	const auto line_vector = *this - other;
	for (size_t i = 0; i < dim; ++i) {
		if (i == axis_tag)
			continue;

		if (std::abs(line_vector[i]) > 1.0E-10)
			return false;
	}
	return true;
}

template <size_t dim> 
std::string EuclideanVector<dim>::to_string(void) const {
	std::string result;
	for (const auto& element : this->small_buffer_)
		result += ms::double_to_string(element) + " ";
	result.pop_back();	
	return result;
}

template <size_t dim> std::ostream& operator<<(std::ostream& os, const EuclideanVector<dim>& x) {
	return os << x.to_string();
};

template <size_t dim>
EuclideanVector<dim> operator*(const double constant, const EuclideanVector<dim>& x) {
	return x * constant;
}