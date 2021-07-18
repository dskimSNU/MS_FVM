#pragma once
#include "Text.h"
#include "Governing_Equation.h"
#include "Spatial_Discrete_Method.h"
#include "Reconstruction_Method.h"
#include "Inital_Condition.h"

#include <iostream>
#include <sstream>


class Log
{
private:
	inline static std::string path_;
	inline static Text log_txt_;

public:
	inline static std::ostringstream content_;

private:
	Log(void) = delete;

public:
	template <typename Governing_Equation, typename Spatial_Discrete_Method, typename Reconstruction_Method, typename Initial_Condition>
	static void initialize(const std::string& grid_file_name) {
		static_require(ms::is_governing_equation<Governing_Equation>, "Wrong governing equation");
		static_require(ms::is_spatial_discrete_method<Spatial_Discrete_Method>, "Wrong spatial discrete method");
		static_require(ms::is_reconsturction_method<Reconstruction_Method>, "Wrong reconstruction method");
		static_require(ms::is_initial_condition<Initial_Condition>, "Wrong initial condition");

		Log::path_ = "Result/Post/" + Governing_Equation::name()
			+ "/" + Spatial_Discrete_Method::name()
			+ "/" + Reconstruction_Method::name()
			+ "/" + Initial_Condition::name()
			+ "/" + grid_file_name + "/_log.txt";
	}

	static void print(void) {
		auto log_str = Log::content_.str();

		std::cout << log_str;
		Log::log_txt_ << std::move(log_str);

		std::ostringstream tmp;
		Log::content_.swap(tmp);
	}

	static void write(void) {
		Log::log_txt_.write(path_);
	}
};

//class Log
//{
//private:
//	std::string path_;
//	std::ostringstream log_stream_;
//	Text log_txt_;	
//
//private:
//	Log(void) = default;
//
//	~Log(void) {
//		log_txt_.write("log.txt");
//	}
//
//public:
//	template <typename T>
//	Log& operator<<(const T& sv) {
//		log_stream_ << sv;
//		return *this;
//	}
//
//	inline void print(void) {
//		static const auto original_flag = std::ostringstream().flags();
//		
//		auto log_str = log_stream_.str();
//		std::cout << log_str;
//		log_txt_ << std::move(log_str);
//
//		log_stream_.clear();
//		log_stream_.flags(original_flag);
//	}
//
//	static Log& instance(void) {
//		static Log instance;
//		return instance;
//	}
//
//	template <typename Governing_Equation, typename Spatial_Discrete_Method, typename Reconstruction_Method, typename Initial_Condition>
//	void initialize(const std::string& grid_file_name) {
//		static_require(ms::is_governing_equation<Governing_Equation>,			"Wrong governing equation");
//		static_require(ms::is_spatial_discrete_method<Spatial_Discrete_Method>, "Wrong spatial discrete method");
//		static_require(ms::is_reconsturction_method<Reconstruction_Method>,		"Wrong reconstruction method");
//		static_require(ms::is_initial_condition<Initial_Condition>,				"Wrong initial condition");
//
//		path_ = "Result/Post/" + Governing_Equation::name()
//			+ "/" + Spatial_Discrete_Method::name()
//			+ "/" + Reconstruction_Method::name()
//			+ "/" + Initial_Condition::name()
//			+ "/" + grid_file_name + "/log.txt";
//	}
//};



