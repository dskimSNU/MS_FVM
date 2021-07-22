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

	Text& operator<<(const std::string& str);
	Text& operator<<(std::string&& str);


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

//#include "FatalError.h"
//#include "FreeFunction.h"

//
//
//class Text
//{
//	using Iter = std::vector<std::string>::iterator;
//	using CIter = std::vector<std::string>::const_iterator;
//
//private:
//	std::vector<std::string> sentence_set_;
//
//public:
//	Text(void) = default;
//	
//	Text(const std::vector<std::string>& other_sentence_set)
//		: sentence_set_(other_sentence_set) {};
//
	//Text(std::initializer_list<std::string> list)
	//	: sentence_set_{ list } {};
//
//	Text(const std::string& read_file);
//
//	Text(std::ifstream& read_file_stream, const size_t num_line = -1);
//
//	Text(std::ifstream& read_file_stream, const std::streampos& start_position, const size_t num_line = -1);
//
//
//
//	template <typename T>
//	Text& operator<<(const T& value) {
//		this->sentence_set_.emplace_back(Editor::to_String(value));
//		return *this;
//	};
//
//	Text& operator<<(std::string&& sentence) {
//		this->sentence_set_.emplace_back(std::move(sentence));
//		return *this;
//	};
//
//	std::string& operator[](const size_t index) {
//		return this->sentence_set_[index];
//	};
//
//	const std::string& operator[](const size_t index) const {
//		return this->sentence_set_[index];
//	};
//
//
//	std::string& back(void) {
//		return this->sentence_set_.back();
//	};
//
//	const std::string& back(void) const {
//		return this->sentence_set_.back();
//	};
//
//	Iter begin(void) {
//		return this->sentence_set_.begin(); 
//	};
//
//	CIter begin(void) const {
//		return this->sentence_set_.begin();
//	};
//
//	Iter end(void) {
//		return this->sentence_set_.end();
//	};
//
//	CIter end(void) const {
//		return this->sentence_set_.end();
//	};
//
//	Iter erase(const Iter& iter) {
//		return this->sentence_set_.erase(iter);
//	};
//
//	std::string& front(void) {
//		return this->sentence_set_.front();
//	};
//
//	const std::string& front(void) const {
//		return this->sentence_set_.front();
//	};
//
//	void pop_back(void) {
//		this->sentence_set_.pop_back();
//	};
//
//	void reserve(const size_t required_capacity) {
//		this->sentence_set_.reserve(required_capacity);
//	};
//
//	size_t size(void) const {
//		return this->sentence_set_.size();
//	};
//
//
//	void add_Write(const std::string& file_path) const;
//
//	void merge(const Text& other) {
//		this->sentence_set_.insert(this->end(), other.begin(), other.end());
//	};
//
//	void merge(Text&& other) {
//		this->sentence_set_.insert(this->end(), std::make_move_iterator(other.begin()), std::make_move_iterator(other.end()));
//	};
//
//	template <typename T>
//	void remove(const T& target) {
//		this->sentence_set_.erase(std::remove(this->begin(), this->sentence_set_.end(), Editor::to_String(target)), this->sentence_set_.end());
//	};
//
//	void remove_empty(void) {
		//this->sentence_set_.erase(std::remove_if(this->begin(), this->end(), [](const std::string& str) {return str.empty(); }), this->end());
//	};
//
//	void replace(const std::string& old_str, const std::string& new_str) {
//		for (std::string& sentence : this->sentence_set_)
//			Editor::replace(sentence, old_str, new_str);
//	};
//
//	void write(const std::string& file_path) const;
//};
//
//

