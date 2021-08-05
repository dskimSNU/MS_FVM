#pragma once
#include <array>
#include <iomanip>
#include <initializer_list>
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
class Euclidean_Vector
{
private:
	std::array<double, dim> vals_ = { 0 };

public:
	Euclidean_Vector(void) = default;
	Euclidean_Vector(const std::array<double, dim>& other) : vals_(other) {};
	template <typename... Args>
	Euclidean_Vector(Args... args);

	Euclidean_Vector& operator+=(const Euclidean_Vector& y);
	Euclidean_Vector& operator-=(const Euclidean_Vector& y);
	Euclidean_Vector& operator*=(const double scalar);
	Euclidean_Vector operator+(const Euclidean_Vector& y) const;	
	Euclidean_Vector operator-(const Euclidean_Vector& y) const;
	Euclidean_Vector operator*(const double scalar) const;
	bool operator==(const Euclidean_Vector& y) const;
	bool operator==(const Euclidean_Vector<0>& y) const;
	double operator[](const size_t position) const;
	double& operator[](const size_t position);

	double at(const size_t position) const;
	Euclidean_Vector& be_normalize(void);
	static constexpr size_t dimension(void);
	const double* data(void) const;
	double inner_product(const Euclidean_Vector& y) const;
	double L1_norm(void) const;
	double norm(void) const;
	bool is_axis_translation(const Euclidean_Vector& other, const size_t axis_tag) const;
	std::string to_string(void) const;
};


//user-defined deduction guides
template <typename... Args>
Euclidean_Vector(Args... args)->Euclidean_Vector<sizeof...(Args)>;

template <size_t dim> 
std::ostream& operator<<(std::ostream& os, const Euclidean_Vector<dim>& x);

template <size_t dim>
Euclidean_Vector<dim> operator*(const double constant, const Euclidean_Vector<dim>& x);


using Dynamic_Euclidean_Vector_ = Euclidean_Vector<0>;

template <>
class Euclidean_Vector<0>
{
private:
	template <size_t dim>
	friend class Euclidean_Vector;

private:
	std::vector<double> vals_;

public:
	Euclidean_Vector(const size_t dimension) : vals_(dimension) {};
	Euclidean_Vector(const std::initializer_list<double> list) : vals_(list) {};

	double& operator[](const size_t position);
	bool operator==(const Dynamic_Euclidean_Vector_& other) const;

	double at(const size_t position) const;
	size_t dimension(void) const;
	std::string to_string(void) const;
};


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
Euclidean_Vector<dim>::Euclidean_Vector(Args... args) : vals_{ static_cast<double>(args)... } {
	static_require(sizeof...(Args) <= dim, "Number of arguments can not exceed dimension");
	static_require(ms::are_arithmetics<Args...>, "every arguments should be arithmetics");
};


template <size_t dim>
Euclidean_Vector<dim>& Euclidean_Vector<dim>::operator+=(const Euclidean_Vector& y) {
	for (size_t i = 0; i < dim; ++i)
		this->vals_[i] += y.vals_[i];
	return *this;
}

template <size_t dim>
Euclidean_Vector<dim>& Euclidean_Vector<dim>::operator-=(const Euclidean_Vector& y) {
	for (size_t i = 0; i < dim; ++i)
		this->vals_[i] -= y.vals_[i];
	return *this;
}

template <size_t dim>
Euclidean_Vector<dim>& Euclidean_Vector<dim>::operator*=(const double scalar) {
	for (size_t i = 0; i < dim; ++i)
		this->vals_[i] *= scalar;
	return *this;
}

template <size_t dim> 
Euclidean_Vector<dim> Euclidean_Vector<dim>::operator+(const Euclidean_Vector& y) const {
	auto result = *this;
	return result += y;
}

template <size_t dim> Euclidean_Vector<dim> Euclidean_Vector<dim>::operator-(const Euclidean_Vector& y) const {
	auto result = *this;
	return result -= y;
}

template <size_t dim> Euclidean_Vector<dim> Euclidean_Vector<dim>::operator*(const double scalar) const {
	auto result = *this;
	return result *= scalar;
}

template <size_t dim> 
bool Euclidean_Vector<dim>::operator==(const Euclidean_Vector& y) const {
	for (size_t i = 0; i < dim; ++i) {
		if (this->vals_[i] != y.vals_[i])
			return false;
	}
	return true;
}

template <size_t dim>
bool Euclidean_Vector<dim>::operator==(const Euclidean_Vector<0>& y) const {
	if (dim != y.dimension())
		return false;
	
	for (size_t i = 0; i < dim; ++i) {
		if (this->vals_[i] != y.vals_[i])
			return false;
	}
	return true;
}

template <size_t dim> 
double Euclidean_Vector<dim>::operator[](const size_t position) const {
	dynamic_require(position <= dim, "Position should be less than dimension");
	return this->vals_[position];
}

template <size_t dim>
double& Euclidean_Vector<dim>::operator[](const size_t position) {
	dynamic_require(position <= dim, "Position should be less than dimension");
	return this->vals_[position];
}


template <size_t dim>
double Euclidean_Vector<dim>::at(const size_t position) const {
	dynamic_require(position <= dim, "Position should be less than dimension");
	return this->vals_[position];
}

template <size_t dim>
Euclidean_Vector<dim>& Euclidean_Vector<dim>::be_normalize(void) {
	const auto scale_factor = 1.0 / this->norm();
	return (*this) *= scale_factor;
}

template <size_t dim>
constexpr size_t Euclidean_Vector<dim>::dimension(void) {
	return dim;
}


template <size_t dim>
const double* Euclidean_Vector<dim>::data(void) const {
	if constexpr (dim != 0)
		return this->vals_.data();
}

template <size_t dim>
double Euclidean_Vector<dim>::inner_product(const Euclidean_Vector& y) const {
	double result = 0;
	for (size_t i = 0; i < dim; ++i)
		result += this->vals_[i] * y.vals_[i];
	return result;
}

template <size_t dim>
double Euclidean_Vector<dim>::L1_norm(void) const {
	double L1_norm = 0.0;
	for (size_t i = 0; i < dim; ++i)
		L1_norm += std::abs(this->vals_[i]);
	return L1_norm;
}

template <size_t dim>
double Euclidean_Vector<dim>::norm(void) const {
	return std::sqrt(this->inner_product(*this));
}

template <size_t dim>
bool Euclidean_Vector<dim>::is_axis_translation(const Euclidean_Vector& other, const size_t axis_tag) const {
	const auto line_vector = *this - other;
	for (size_t i = 0; i < dim; ++i) {
		if (i == axis_tag)
			continue;

		if (std::abs(line_vector.at(i) > 1.0E-10))
			return false;
	}
	return true;
}

template <size_t dim> 
std::string Euclidean_Vector<dim>::to_string(void) const {
	std::string result;
	for (const auto& element : this->vals_)
		result += ms::double_to_string(element) + " ";
	result.pop_back();	
	return result;
}

template <size_t dim> std::ostream& operator<<(std::ostream& os, const Euclidean_Vector<dim>& x) {
	return os << x.to_string();
};

template <size_t dim>
Euclidean_Vector<dim> operator*(const double constant, const Euclidean_Vector<dim>& x) {
	return x * constant;
}