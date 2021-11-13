#pragma once
#include <string>
#include <vector>
#include <unordered_map>
namespace misc
{
    void trim_leading_and_ending_whitespace(std::string& string);
    std::vector<std::string> split_on(const std::string& in, char);
    bool is_string_in_string_vector(std::string, std::vector<std::string>);
    std::string to_lower(std::string);
    bool compare_strings(std::string);
    bool is_keyword(std::string);
    bool is_valid_name(std::string);
    bool is_whitespace(std::string);
    bool is_comment(const std::string&);
    void append_newline_to_textfile(std::string filename);
    std::unordered_map<std::string,std::string> options_map(std::string line);
}