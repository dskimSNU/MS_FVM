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
template <size_t dimension>
class EuclideanVector
{
public:
	EuclideanVector(void) = default;
	EuclideanVector(const std::array<double, dimension>& other) : small_buffer_(other) {};
	EuclideanVector(const std::vector<double>& other) : elements_(other) {};
	template <typename... Args>
	EuclideanVector(Args... args);

	EuclideanVector& operator+=(const EuclideanVector& y);
	EuclideanVector& operator-=(const EuclideanVector& y);
	EuclideanVector& operator*=(const double scalar);
	EuclideanVector operator+(const EuclideanVector& y) const;	
	EuclideanVector operator-(const EuclideanVector& y) const;
	EuclideanVector operator*(const double scalar) const;
	bool operator<(const EuclideanVector& y) const;
	bool operator==(const EuclideanVector& y) const;
	double operator[](const size_t position) const;

	EuclideanVector& be_normalize(void);
	//static constexpr size_t dimension(void);
	const double* data(void) const;
	double inner_product(const EuclideanVector& y) const;
	double norm(void) const;	
	std::string to_string(void) const;

private:
	std::array<double, dimension> small_buffer_ = { 0 };
	std::vector<double> elements_;
};


//user-defined deduction guides
template <typename... Args>
EuclideanVector(Args... args)->EuclideanVector<sizeof...(Args)>;
EuclideanVector(const std::vector<double>& vec)->EuclideanVector<0>;


template <size_t dimension> 
std::ostream& operator<<(std::ostream& os, const EuclideanVector<dimension>& x);

template <size_t dimension>
EuclideanVector<dimension> operator*(const double constant, const EuclideanVector<dimension>& x);


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
template <size_t dimension>
template <typename... Args>
EuclideanVector<dimension>::EuclideanVector(Args... args) : small_buffer_{ static_cast<double>(args)... } {
	static_require(sizeof...(Args) <= dimension, "Number of arguments can not exceed dimension");
	static_require(ms::are_arithmetics<Args...>, "every arguments should be arithmetics");
};


template <size_t dimension>
EuclideanVector<dimension>& EuclideanVector<dimension>::operator+=(const EuclideanVector& y) {
	for (size_t i = 0; i < dimension; ++i)
		this->small_buffer_[i] += y.small_buffer_[i];
	return *this;
}

template <size_t dimension>
EuclideanVector<dimension>& EuclideanVector<dimension>::operator-=(const EuclideanVector& y) {
	for (size_t i = 0; i < dimension; ++i)
		this->small_buffer_[i] -= y.small_buffer_[i];
	return *this;
}

template <size_t dimension>
EuclideanVector<dimension>& EuclideanVector<dimension>::operator*=(const double scalar) {
	for (size_t i = 0; i < dimension; ++i)
		this->small_buffer_[i] *= scalar;
	return *this;
}

template <size_t dimension> 
EuclideanVector<dimension> EuclideanVector<dimension>::operator+(const EuclideanVector& y) const {
	auto result = *this;
	return result += y;
}

template <size_t dimension> EuclideanVector<dimension> EuclideanVector<dimension>::operator-(const EuclideanVector& y) const {
	auto result = *this;
	return result -= y;
}

template <size_t dimension> EuclideanVector<dimension> EuclideanVector<dimension>::operator*(const double scalar) const {
	auto result = *this;
	return result *= scalar;
}

template <size_t dimension>
bool EuclideanVector<dimension>::operator<(const EuclideanVector& y) const {
	for (size_t i = 0; i < dimension; ++i) {
		if (this->small_buffer_[i] == y.small_buffer_[i])
			continue;
		return this->small_buffer_[i] < y.small_buffer_[i];
	}
	return false;
}

template <size_t dimension> 
bool EuclideanVector<dimension>::operator==(const EuclideanVector& y) const {
	for (size_t i = 0; i < dimension; ++i) {
		if (this->small_buffer_[i] != y.small_buffer_[i])
			return false;
	}
	return true;
}

template <size_t dimension> 
double EuclideanVector<dimension>::operator[](const size_t position) const {
	dynamic_require(position <= dimension, "Position can not exceed dimension");
	return this->small_buffer_[position];
}

template <size_t dimension>
EuclideanVector<dimension>& EuclideanVector<dimension>::be_normalize(void) {
	const auto scale_factor = 1.0 / this->norm();
	return (*this) *= scale_factor;
}

//template <size_t dimension>
//constexpr size_t EuclideanVector<dimension>::dimension(void) {
//	return dimension;
//}


template <size_t dimension>
const double* EuclideanVector<dimension>::data(void) const {
	if constexpr (dimension != 0)
		return this->small_buffer_.data();
	else
		return this->elements_.data();
}

template <size_t dimension>
double EuclideanVector<dimension>::inner_product(const EuclideanVector& y) const {
	double result = 0;
	for (size_t i = 0; i < dimension; ++i)
		result += this->small_buffer_[i] * y.small_buffer_[i];
	return result;
}

template <size_t dimension>
double EuclideanVector<dimension>::norm(void) const {
	return std::sqrt(this->inner_product(*this));
}

template <size_t dimension> 
std::string EuclideanVector<dimension>::to_string(void) const {
	std::string result;
	for (const auto& element : this->small_buffer_)
		result += ms::double_to_string(element) + " ";
	result.pop_back();	
	return result;
}

template <size_t dimension> std::ostream& operator<<(std::ostream& os, const EuclideanVector<dimension>& x) {
	return os << x.to_string();
};

template <size_t dimension>
EuclideanVector<dimension> operator*(const double constant, const EuclideanVector<dimension>& x) {
	return x * constant;
}