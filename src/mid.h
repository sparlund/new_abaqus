#pragma once
#include <string>
#include <vector>
#include <unordered_map>
		
class Mid
{
private:
    /* data */
public:
    // std::array<float, 3> color;
    unsigned int id;
    std::string name;
    
    Mid(int id, std::string name);
    std::unordered_map<std::string,std::string> options;
    ~Mid();
};
