#pragma once
#include <vector>
#include <string>
#include <memory>
#include <iostream>

template <class T>
class Set
{
private:
    std::string name;
    std::vector<std::shared_ptr<T>> entities;
public:
    unsigned int get_number_of_entities(){return entities.size();}
    std::string get_set_name(){return name;}
    std::shared_ptr<T> get_entity(unsigned int i){return entities.at(i);}
    void add_entity(std::shared_ptr<T> entity_pointer){entities.push_back(entity_pointer);}
    Set(std::string name):name(name){
        std::cout << "*NSET: nset = " << name;
    };
    ~Set(){};
};





