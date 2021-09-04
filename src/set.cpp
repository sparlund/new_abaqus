#include "element.h"
#include "node.h"
#include "misc_string_functions.h"
#include "set.h"

template class Set<Node*>;
template class Set<Element*>;

template <class T>
void Set<T>::add_entities(std::string line){
    // split on comma
    // std::vector<std::string> entities = misc::split_on(std::string in, ',');
    // for (const auto& entity : entities)
    // {

    // }
};

template <class T>
void Set<T>::add_entity(T entity_pointer){
    entities.emplace_back(entity_pointer);
    }
template <class T>
Set<T>::Set(const std::string& name):name{name}{
        std::cout << "*NSET: nset = " << name;
    };
