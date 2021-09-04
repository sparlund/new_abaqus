#pragma once
#include <iostream>
#include <memory>
#include <string>
#include <vector>


// NOTE: T is pointer!
template <class T>
class Set
{
private:
    const std::string name;
    std::vector<T> entities;
public:
    size_t get_number_of_entities() const {return entities.size();}
    std::string get_set_name() const {return name;}
    T get_entity(size_t i){return entities.at(i);}
    void add_entities(std::string line);
    void add_entity(T entity_pointer);
    Set(std::string name);
    ~Set() = default;
};