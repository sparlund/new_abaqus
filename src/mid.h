#pragma once
#include <string>
#include <vector>
#include <unordered_map>
		
class Mid
{
private:
    /* data */
public:
    unsigned int id;
    std::string name;
    float density;
    Mid(int id, std::string name);
    std::unordered_map<std::string,std::string> options;
    ~Mid();
};
