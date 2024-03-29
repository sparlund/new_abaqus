#include "../include/set.h"
#include "../include/element.h"
#include "../include/misc_string_functions.h"
#include "../include/node.h"

template class Set<Node*>;
template class Set<Element*>;

template <class T>
void Set<T>::add_entity(T entity){
    entities.emplace_back(entity);
    }
template <class T>
Set<T>::Set(const std::string& name):name{name}{
        std::cout << "*NSET: nset = " << name << std::endl;
    };
template <class T>
bool Set<T>::is_entity_in_set(T in) const
{
    for(const auto& entity: entities)
    {
        if(in->id == entity->id)
        {
            return true;
        }
    }
    return false;
}
