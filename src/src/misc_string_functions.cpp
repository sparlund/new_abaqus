#include <algorithm>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include "../include/misc_string_functions.h"
#include <unordered_map>

bool misc::is_string_in_string_vector(std::string input_string, std::vector<std::string> input_vector){
    for (unsigned int i = 0; i < input_vector.size(); i++)
    {
        if (input_string == input_vector.at(i))
        {
            return true;
        }
        
    }
    return false;
    
}

bool misc::is_valid_name(std::string entity_name){
    if (entity_name.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ01234567890_;") != std::string::npos)
    {
        return false;
    }
    else{
        return true;
    }
};

void misc::trim_leading_and_ending_whitespace(std::string& string){
    const auto strBegin = string.find_first_not_of(" ");
    if (strBegin == std::string::npos){
        return;
    }
    const auto strEnd = string.find_last_not_of(" ");
    const auto strRange = strEnd - strBegin + 1;
    string = string.substr(strBegin, strRange);
    return;
};

std::unordered_map<std::string,std::string> misc::options_map(std::string line){
    // returns a map with options as keys and the option values as parameters
    std::unordered_map<std::string, std::string> options;
    // split input line on comma, keep as capital case
    auto strings = misc::split_on(line,',');
    // loop over comma-split string and start creating 
    // key-value pairs
    for (unsigned int i = 0; i < strings.size(); i++)
    {
        // try to split on equal sign
        auto key_and_value = misc::split_on(strings.at(i), '=');
        auto key = key_and_value.at(0);
        // trim leading and trailing whitespace of string
        misc::trim_leading_and_ending_whitespace(key);
        if(key_and_value.size() == 1)
        {
            options[key] = key;
        }
        else
        {
            auto value = key_and_value.at(1);
            options[key] = value;
        }        
    }
    return options;
};

void misc::append_newline_to_textfile(std::string filename)
{
    std::ofstream input_file(filename, std::ios::app);
    
}

bool misc::is_keyword(std::string line){
    // If the line starts with an asterix
    // and the second character is an alpha
    // we don't have a keyword!
    // NOTE: compare return 0 if equal, isalpha returns non-zero if true
    if (line.compare(0, 1, "*") == 0 && isalpha(line.at(1))!= 0){
        return true;
    }
    else{    
        return false;
    }
};

bool misc::is_comment(const std::string& line){
    // NOTE: compare return 0 if equal
    if (line.compare(0, 2, "**") == 0){
        return true;
    }
    else {
    return false;
    };   

};

bool misc::is_whitespace(std::string s){
    for(size_t index = 0; index < s.length(); index++){
        if(!std::isspace(s[index])){
            return false;
        }
    }
    return true;
}

std::vector<std::string> misc::split_on(const std::string& in,char split_char)
{
    std::stringstream ss(in);
    std::vector<std::string> items;
    while (ss.good())
    {
        std::string item;
        getline(ss, item, split_char);
        if(!is_whitespace(item)){
            items.push_back(item);
        }
    };
    return items;
};



std::string misc::to_lower(std::string text)
{
    // return a string in all lower case
    std::for_each(text.begin(), text.end(), [](char &c) {
        c = std::tolower(c);
    });
    return text;
};