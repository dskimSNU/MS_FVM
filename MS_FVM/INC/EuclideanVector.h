#pragma once
#include <array>
#include <iomanip>
#include <string>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <vector>

#define dynamic_require(requirement, state) if (!(requirement)) throw std::runtime_error(state)

namespace ms {
	template <size_t size1, size_t size2> using is_same = std::enable_if_t<(size1 == size2), bool>;
	template <typename... Args> using are_arithmetics = std::enable_if_t<std::conjunction_v<std::is_arithmetic<Args>...>, bool>;
}

// CLASS TEMPLATE euclidean_vector
template <size_t small = 0>
class EuclideanVector
{
public:
	EuclideanVector(void) = default;
	template <typename... Args, ms::is_same<sizeof...(Args), small> = true, ms::are_arithmetics<Args...> = true>
	EuclideanVector(Args... args) : small_buffer_{ static_cast<double>(args)... } {};
	EuclideanVector(const std::array<double, small>& other) : small_buffer_(other) {};
	EuclideanVector(const std::vector<double>& other) : elements_(other) {};

	EuclideanVector operator+(const EuclideanVector& y) const;
	EuclideanVector operator-(const EuclideanVector& y) const;
	EuclideanVector operator*(const double scalar) const;
	bool operator==(const EuclideanVector& y) const;
	double operator[](const size_t position) const;

	//constexpr size_t dimension(void) const;
	double inner_product(const EuclideanVector& y) const;
	double norm(void) const;
	std::string to_string(void) const;
	

private:
	std::array<double, small> small_buffer_ = { 0 };
	std::vector<double> elements_;	//template specialization when small =0
};

template <typename... Args> 
EuclideanVector(Args... args)->EuclideanVector<sizeof...(Args)>;  //user-defined deduction guides

template <size_t small> 
std::ostream& operator<<(std::ostream& os, const EuclideanVector<small>& x);

template <size_t small>
EuclideanVector<small> operator*(const double constant, const EuclideanVector<small>& x);


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
template <size_t small> EuclideanVector<small> EuclideanVector<small>::operator+(const EuclideanVector& y) const {
	auto result = *this;
	for (size_t i = 0; i < small; ++i)
		result.small_buffer_[i] += y.small_buffer_[i];
	return result;
}

template <size_t small> EuclideanVector<small> EuclideanVector<small>::operator-(const EuclideanVector& y) const {
	auto result = *this;
	for (size_t i = 0; i < small; ++i)
		result.small_buffer_[i] -= y.small_buffer_[i];
	return result;
}

template <size_t small> EuclideanVector<small> EuclideanVector<small>::operator*(const double scalar) const {
	auto result = *this;
	for (auto& element : result.small_buffer_)
		element *= scalar;
	return result;
}

template <size_t small> bool EuclideanVector<small>::operator==(const EuclideanVector& y) const {
	for (size_t i = 0; i < small; ++i) {
		if (this->small_buffer_[i] != y.small_buffer_[i])
			return false;
	}
	return true;
}

template <size_t small> double EuclideanVector<small>::operator[](const size_t position) const {
	dynamic_require(position <= small, "EuclideanVector out of range");
	return this->small_buffer_[position];
}

template <size_t small>
double EuclideanVector<small>::inner_product(const EuclideanVector& y) const {
	double result = 0;
	for (size_t i = 0; i < small; ++i)
		result += this->small_buffer_[i] * y.small_buffer_[i];
	return result;
}

template <size_t small>
double EuclideanVector<small>::norm(void) const {
	return std::sqrt(this->inner_product(*this));
}

template <size_t small> 
std::string EuclideanVector<small>::to_string(void) const {
	std::string result;
	for (const auto& element : this->small_buffer_)
		result += ms::double_to_string(element) + " ";
	result.pop_back();	
	return result;
}

template <size_t small> std::ostream& operator<<(std::ostream& os, const EuclideanVector<small>& x) {
	return os << x.to_string();
};

template <size_t small>
EuclideanVector<small> operator*(const double constant, const EuclideanVector<small>& x) {
	return x * constant;
}














//template <size_t small> EuclideanVector<small> EuclideanVector<small>::operator+(const EuclideanVector& y) const {
//	return this->vector_ + y.vector_;
//}
//template <size_t small> EuclideanVector<small> EuclideanVector<small>::operator-(const EuclideanVector& y) const {
//	return this->vector_ - y.vector_;
//}
//template <size_t small> EuclideanVector<small> EuclideanVector<small>::operator*(const double scalar) const {
//	return this->vector_ * scalar;
//}
//template <size_t small> bool EuclideanVector<small>::operator==(const EuclideanVector& y) const {
//	return this->vector_ == y.vector_;
//}
//template <size_t small> constexpr size_t EuclideanVector<small>::dimension(void) const {
//	return small;
//}
//template <size_t small> std::string EuclideanVector<small>::to_string(void) const {
//	return this->vector_.to_string();
//}
//
//
//template <size_t small>
//std::ostream& operator<<(std::ostream& os, const EuclideanVector<small>& x) {
//	return os << x.to_string();
//};

//#include <algorithm>
//#include <fstream>
//#include <iomanip>	//set precision
//#include <mkl.h>
//#include <stdexcept>
//#include <string>
//#include <sstream>
//#include <vector>
//
////variable vector //real vector
//class MathVector : public std::vector<double>
//{
//public:
//	template <typename ... Vals>
//	explicit MathVector(Vals&&... values) :std::vector<double>(std::forward<Vals>(values)...) {};
//	MathVector(std::initializer_list<double> list) :std::vector<double>(list) {};
//	MathVector(const std::vector<double>& vec) :std::vector<double>(vec) {};
//	MathVector(const double) = delete;
//	
//	MathVector& operator+=(const MathVector& y);
//	MathVector& operator-=(const MathVector& y);
//	MathVector& operator*=(const double scalar);
	//MathVector operator+(const MathVector& y) const;
	//MathVector operator-(const MathVector& y) const;
	//MathVector operator*(const double scalar) const;
//
//	MathVector& abs(void);	
//	bool compare_with_finite_precision(const MathVector& other, const size_t ULP_precision = 4) const;
//	double inner_product(const MathVector& other) const;
//	double L2_Norm(void) const;
//	void merge(const MathVector& other);
//	void merge(MathVector&& other);
//	MathVector& normalize(void);
//	std::string to_string(void) const;
//};
//std::ostream& operator<<(std::ostream& os, const MathVector& x);
//MathVector operator*(const double scalar, const MathVector& x);
//
//
//namespace ms {
//	MathVector abs(const MathVector& x);
//	MathVector normalize(const MathVector& x);
//	std::string double_to_string(const double val, const size_t precision = 17);
//	bool compare_double(const double d1, const double d2, const size_t ULP_precision = 4);
//}
//
//
//
//
//
//template <typename T>
//class VectorFunction : public std::vector<T> // template 특수화를 통해 하나로 묶어보자
//{
//public:
//	template <typename ... Vals>
//	explicit VectorFunction(Vals&&... values) :std::vector<T>(std::forward<Vals>(values)...) {};
//	VectorFunction(std::initializer_list<T> list) :std::vector<T>(list) {};
//
//	MathVector operator()(const MathVector& variable_vector) const {
//		MathVector result;
//		result.reserve(this->size());
//		for (const auto& function : *this)
//			result.push_back(function(variable_vector));
//		return result;
//	}
//
//	std::vector<MathVector> operator()(const std::vector<MathVector>& variable_vector_set) const {
//		std::vector<MathVector> result;
//		result.reserve(variable_vector_set.size());
//		for (const auto& variable_vector : variable_vector_set)
//			result.push_back((*this)(variable_vector));
//		return result;
//	}
//
//	VectorFunction<T> cross_product(const VectorFunction<T>& other) const {
//		constexpr size_t range_dimension = 3;
//		if (this->size() != range_dimension || other.size() != range_dimension)
//			throw std::runtime_error("cross product is only defined on R^3 range");
//
//		VectorFunction<T> result(range_dimension, 0.0);
//		result[0] = (*this)[1] * other[2] - (*this)[2] * other[1];
//		result[1] = (*this)[2] * other[0] - (*this)[0] * other[2];
//		result[2] = (*this)[0] * other[1] - (*this)[1] * other[0];
//
//		return result;
//	}
//
//	VectorFunction<T>& be_derivative(const size_t variable_index) {
//		for (auto& func : *this)
//			func.be_derivative(variable_index);
//		return *this;
//	}
//
//	VectorFunction<T> differentiate(const size_t variable_index) const {
//		auto result = *this;
//		return result.be_derivative(variable_index);
//	}
//
//	size_t domain_dimension(void) const {
//		const size_t num_function = this->size();
//
//		std::vector<size_t> domain_dimension_set(num_function);
//		for (size_t i = 0; i < num_function; ++i)
//			domain_dimension_set[i] = (*this)[i].domain_dimension();
//
//		return *std::max_element(domain_dimension_set.begin(), domain_dimension_set.end());
//	}
//
//	//T L2_norm(void) const {
//	//	T result;
//	//	for (const auto& func : *this)
//	//		result += (func ^ 2);
//	//	return result.power(0.5);
//	//}
//
//	std::string to_string(void) const {
//		std::string result = "{ ";
//		for (const auto& func : *this)
//			result += func.to_string() + ", ";
//		result.erase(result.end() - 2, result.end());
//		result += " }";
//		return result;
//	}
//};

//template <typename T>
//std::ostream& operator<<(std::ostream& os, const VectorFunction<T>& x) {
//	return os << x.to_string();
//};

//template <size_t small>
//class EuclideanVector
//{
//public:
//	template <typename... Args> EuclideanVector(Args... args) : vector_(static_cast<double>(args)...) {};
//	EuclideanVector(const Vector<double, small>& vector) : vector_(vector) {};
//
//	EuclideanVector operator+(const EuclideanVector& y) const;
//	EuclideanVector operator-(const EuclideanVector& y) const;
//	EuclideanVector operator*(const double scalar) const;
//	bool operator==(const EuclideanVector& y) const;
//		
//	constexpr size_t dimension(void) const;
//	std::string to_string(void) const;
//
//private:
//	Vector<double, small> vector_;
//};
//template <typename... Args> EuclideanVector(Args... args)->EuclideanVector<sizeof...(Args)>; //user-defined deduction guides
//template <size_t small> std::ostream& operator<<(std::ostream& os, const EuclideanVector<small>& x);

//template <typename T, size_t small>
//class Vector
//{
//protected:
//	template <typename... Args>
//	Vector(Args... args) : elements_{ args... } {};
//
//public:
//	Vector operator+(const Vector& y) const;
//	Vector operator-(const Vector& y) const;
//	Vector operator*(const double scalar) const;
//	bool operator==(const Vector& y) const;
//
//	constexpr size_t dimension(void) const;
//	std::string to_string(void) const;
//
//protected:
//	std::array<T, small> elements_;
//};
//template <typename T, size_t small>
//std::ostream& operator<<(std::ostream& os, const Vector<T, small>& x);
//
//
//template <size_t small>
//class EuclideanVector : public Vector<double, small>
//{
//public:
//	template <typename... Args>
//	EuclideanVector(Args... args) : Vector<double, sizeof...(Args)>(static_cast<double>(args)...) {};
//
//};
//template <typename... Args> //user-defined deduction guides
//EuclideanVector(Args... args)->EuclideanVector<sizeof...(Args)>; 
//template <size_t small>
//std::ostream& operator<<(std::ostream& os, const EuclideanVector<small>& x);
//
//namespace ms {
//	template <typename T>
//	std::string to_string(const T& arg);
//}