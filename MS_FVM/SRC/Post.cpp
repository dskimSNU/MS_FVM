#include "../INC/Post.h"

Text Post::header_text(const Post_File_Type file_type, const double current_time) {
	static size_t strand_id = 0;

	Text header;
	header.reserve(10);
	header << "Title = " + path_;
	if (file_type == Post_File_Type::Grid) {
		header << "FileType = Grid";
		header << grid_variables_;
	}
	else {
		header << "FileType = Solution";
		header << solution_variables_;
		strand_id++;
	}
	header << "Zone T = " + path_;
	header << zone_type_;
	header << "Nodes = " + std::to_string(num_node_);
	header << "Elements = " + std::to_string(num_element_);
	header << "DataPacking = Block";
	header << "StrandID = " + std::to_string(strand_id);
	header << "SolutionTime = " + ms::double_to_string(current_time) + "\n";

	return header;
}