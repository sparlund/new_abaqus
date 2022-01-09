#pragma once
#include <iostream>
#include <memory>
#include <string>
#include <vector>


// NOTE: T is pointer!
template <class T>
class Set
{
public:
    const std::string name;
    std::vector<T> entities;
    bool is_entity_in_set(T) const;
    size_t size() const {return entities.size();}
    T get_entity(size_t i) const {return entities.at(i);}
    void add_entity(T entity_pointer);
    Set<T>(const std::string& name);
    Set<T>(Set<T>&) = default;
};