#include "set.h"
#include "elements/element.h"
#include "misc_string_functions.h"
#include "node.h"

template class Set<Node*>;
template class Set<Element*>;

template <class T>
void Set<T>::add_entity(T entity_pointer){
    entities.emplace_back(entity_pointer);
    }
template <class T>
Set<T>::Set(const std::string& name):name{name}{
        std::cout << "*NSET: nset = " << name << std::endl;
    };
