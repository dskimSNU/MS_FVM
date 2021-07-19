#pragma once
#include "Text.h"

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
	static void set_path(const std::string& path) {	
		Log::path_ = path;
	}

	static void print(void) {
		auto log_str = Log::content_.str();

		std::cout << log_str;
		Log::log_txt_ << std::move(log_str);

		std::ostringstream tmp;
		Log::content_.swap(tmp);
	}

	static void write(void) {
		const auto file_path = Log::path_ + "_log.txt";
		Log::log_txt_.write(file_path);
	}
};