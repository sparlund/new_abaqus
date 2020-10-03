#pragma once
#include <vector>
#include <string>
#include <unordered_map>
namespace misc
{
void trim(std::string& string);
std::vector<std::string> split_on(std::string in, char);
std::string to_lower(std::string);
bool compare_strings(std::string);
bool is_keyword(std::string);
bool is_comment(std::string);
void append_newline_to_textfile(std::string filename);
std::unordered_map<std::string,std::string> options_map(std::string line);
}