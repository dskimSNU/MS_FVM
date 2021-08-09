#pragma once
#include <algorithm>
#include <fstream>	//file stream
#include <filesystem>
#include <stdexcept>
#include <string>
#include <string_view>
#include <sstream>
#include <vector>

#define dynamic_require(requirement, state) if (!(requirement)) throw std::runtime_error(state)

class Text : public std::vector<std::string>
{
public:
	template <typename ... Vals>
	explicit Text(Vals&&... values) : std::vector<std::string>(std::forward<Vals>(values)...) {};
	Text(std::initializer_list<std::string> list) : std::vector<std::string>( list ) {};
	Text(std::ifstream& file, const size_t num_read_line);

	Text& operator<<(const std::string& str);	//Lvalue ref
	Text& operator<<(std::string&& str);		//Rvalue ref


	void add_write(const std::string& write_file_path) const;
	Text& read_line_by_line(const std::string& read_file_path);
	void read(std::ifstream& file, const size_t num_read_line);
	Text& remove_empty_line(void);
	void write(const std::string& write_file_path) const;

private:
	void make_path(std::string_view file_path) const;
};

std::ostream& operator<<(std::ostream& ostream, const Text& text);


namespace ms {
	std::vector<std::string> parse(const std::string& str, const char delimiter);
	template<typename T>
	T string_to_value(const std::string& str);
	template<typename T>
	std::vector<T> string_to_value_set(const std::vector<std::string>& str_set);
	std::string erase(const std::string& str, const std::string& target);
	std::string upper_case(const std::string& str);
	size_t find_icase(const std::string& str, const std::string& target);
	bool is_there_icase(const std::string& str, const std::string& target);
	std::string double_to_str_sp(const double value); //double to string with show point
	Text extract_file_name_text(const std::string& path);
}


//template definition
namespace ms {
	template<typename T>
	T string_to_value(const std::string& str) {
		std::istringstream iss(str);
		T value;
		iss >> value;
		return value;
	};
	template<typename T>
	std::vector<T> string_to_value_set(const std::vector<std::string>& str_set) {
		std::vector<T> result;
		result.reserve(str_set.size());

		for (const auto& str : str_set)
			result.push_back(ms::string_to_value<T>(str));

		return result;
	};
}



