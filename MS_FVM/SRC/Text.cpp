#include "../INC/Text.h"

Text::Text(std::ifstream& file, const size_t num_read_line) {
	this->read(file, num_read_line);
}

Text& Text::operator<<(const std::string& str) {
	this->push_back(str);
	return *this;
}

Text& Text::operator<<(std::string&& str) {
	this->push_back(std::move(str));
	return *this;
}

void Text::add_write(const std::string& file_path) const {
	this->make_path(file_path);
	std::ofstream output_file(file_path, std::ios::app);
	dynamic_require(output_file.is_open(), "Fail to open file" + file_path);

	const auto num_sentence = this->size();
	for (auto i = this->begin(); i != this->end() - 1; ++i)
		output_file << *i << "\n";
	output_file << this->back();

	output_file.close();
}

Text& Text::read_line_by_line(const std::string& file_path) {	
	std::ifstream file_stream(file_path);
	dynamic_require(file_stream.is_open(), "Fail to open file" + file_path);

	std::string str;
	while (std::getline(file_stream, str))
		this->push_back(std::move(str));

	file_stream.close();
	return *this;
}

void Text::read(std::ifstream& file_stream, const size_t num_read_line) {
	dynamic_require(file_stream.is_open(), "Fail to open file");

	size_t index = 0;
	std::string str;
	while (std::getline(file_stream, str)) {
		this->push_back(std::move(str));
		if (++index == num_read_line)
			break;
	}

	//file_stream.close();
}

Text& Text::remove_empty_line(void) {
	this->erase(std::remove(this->begin(), this->end(), ""), this->end());
	return *this;
}

void Text::write(const std::string& file_path) const {
	this->make_path(file_path);
	std::ofstream output_file(file_path);
	dynamic_require(output_file.is_open(), "Fail to open file" + file_path); 

	const auto num_sentence = this->size();
	for (auto i = this->begin(); i != this->end() - 1; ++i)
		output_file << *i << "\n";
	output_file << this->back();

	output_file.close();
}

void Text::make_path(std::string_view file_path) const {
	const auto file_name_size = file_path.size() - file_path.find_last_of("/") - 1;
	file_path.remove_suffix(file_name_size);
	
	if (file_path.empty())
		return;

	std::filesystem::path p(file_path);
	if (std::filesystem::exists(p))
		return;
	else
		std::filesystem::create_directories(p);	
}

namespace ms {
	std::vector<std::string> parse(const std::string& str, const char delimiter) {
		if (str.empty())
			return std::vector<std::string>();

		std::vector<std::string> parsed_string_set;
		for (auto iter1 = str.begin();;) {
			auto iter2 = std::find(iter1, str.end(), delimiter);

			if (iter1 != iter2)
				parsed_string_set.emplace_back(iter1, iter2);

			if (iter2 == str.end())
				return parsed_string_set;

			iter1 = iter2 + 1;
		}
	};

	std::string erase(const std::string& str, const std::string& target) {
		const auto target_size = target.size();

		auto result = str;
		while (true) {
			const auto position = result.find(target);

			if (position == std::string::npos)
				return result;

			result.erase(position, target_size);
		}
	}

	size_t find_icase(const std::string& str, const std::string& target) {
		auto u_str = ms::upper_case(str);
		auto u_target = ms::upper_case(target);

		return u_str.find(u_target);
	}

	std::string upper_case(const std::string& str) {
		auto result = str;
		std::transform(result.begin(), result.end(), result.begin(), toupper);
		//std::transform(result.begin(), result.end(), result.begin(), std::toupper); //filesystem이랑 있으면 충돌
		return result;
	}

	bool is_there_icase(const std::string& str, const std::string& target) {
		return ms::find_icase(str, target) != std::string::npos;
	}

	std::string double_to_str_sp(const double value) {
		std::ostringstream os;
		os << std::setprecision(16) << std::showpoint << value;
		return os.str();
	}
}


std::ostream& operator<<(std::ostream& ostream, const Text& text) {
	for (const auto& line : text)
		ostream << line << "\n";
	return ostream;
}

//	std::vector<std::string> parse(const std::string& str, const std::vector<char>& delimiter_set) {
//		if (str.empty())
//			return std::vector<std::string>();
//		if (delimiter_set.empty())
//			return std::vector<std::string>(1, str);
//
//		auto temporary_string = str;
//
//		const auto& reference_delimiter = delimiter_set.front();
//		for (const auto& other_delimiter : delimiter_set)
//		{
//			if (other_delimiter != reference_delimiter)
//				Editor::replace(temporary_string, other_delimiter, reference_delimiter);
//		}
//
//		return StringEditor::parse(temporary_string, reference_delimiter);
//	};
//
//	std::string& remove_comment(std::string& str, const std::string& comment) {
//		const auto position = Tool::find_First_Position(str, comment);
//
//		if (position == std::string::npos)
//			return str;
//		else
//			return str.erase(position, std::string::npos);
//	}
//
//	std::string& UpperCase(std::string& str) {
//		std::transform(str.begin(), str.end(), str.begin(), std::toupper);
//		return str;
//	};
//
//	std::string UpperCase(const std::string& str) {
//		auto result = str;
//		return StringEditor::UpperCase(result);
//	};
//
//	std::string UpperCase(std::string&& str) {
//		auto result = std::move(str);
//		return StringEditor::UpperCase(result);
//	};
//}
//
//
//Text::Text(std::ifstream& read_file_stream, const size_t num_line) {
//	if (!read_file_stream.is_open() || num_line == 0)
//		return;
//	
//	std::string sentence;
//
//	size_t line_count = 1;
//	while (std::getline(read_file_stream, sentence)) {
//		this->sentence_set_.emplace_back(std::move(sentence));
//		
//		if (line_count++ == num_line)
//			break;
//	}
//}
//
//Text::Text(std::ifstream& read_file_stream, const std::streampos& start_position, const size_t num_line) {
//	if (!read_file_stream.is_open() || num_line == 0)
//		return;
//
//	read_file_stream.seekg(start_position);
//
//	std::string sentence;
//
//	size_t line_count = 1;
//	while (std::getline(read_file_stream, sentence)) {
//		this->sentence_set_.emplace_back(std::move(sentence));
//
//		if (line_count++ == num_line)
//			break;
//	}
//}
//
//void Text::add_Write(const std::string& file_path) const{
	//std::ofstream outfile(file_path, std::ios::app);

	//if (!outfile.is_open())
	//	FATAL_ERROR("Fail to open" + file_path);

	//const auto num_sentence = this->sentence_set_.size();
	//for (size_t i = 0; i < num_sentence - 1; ++i)
	//	outfile << this->sentence_set_[i] << "\n";
	//outfile << this->sentence_set_[num_sentence - 1];

//	outfile.close();
//}
//
//void Text::write(const std::string& file_path) const{
	//std::ofstream outfile(file_path);

	//if (!outfile.is_open())
	//	FATAL_ERROR("Fail to open" + file_path);

	//const auto num_sentence = this->sentence_set_.size();
	//for (size_t i = 0; i < num_sentence - 1; ++i)
	//	outfile << this->sentence_set_[i] << "\n";
	//outfile << this->sentence_set_[num_sentence - 1];

	//outfile.close();
//}