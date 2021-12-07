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
    bool is_entity_in_set(T) const;
    size_t size() const {return entities.size();}
    std::string get_set_name() const {return name;}
    T get_entity(size_t i){return entities.at(i);}
    void add_entity(T entity_pointer);
    Set<T>(const std::string& name);
    Set<T>(Set<T>&) = default;
};