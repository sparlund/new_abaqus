#include "../include/contact.h"
#include "../include/element.h"

Contact::Contact(Set<Node*>& master, Set<Node*>& slave) : master(master), slave(slave)
{
    // Create segments from nodes in master set
    for(size_t i = 0; i < master.get_number_of_entities(); i++)
    {
        auto node = master.get_entity(i);
        auto connected_elements = node->connected_elements;
        for(auto element : connected_elements)
        {
            auto segments = element->get_segments(node);
            for(const auto segment : segments)
            {
                if(master.is_entity_in_set(segment.first) && master.is_entity_in_set(segment.second))
                {
                    master_segments.push_back(segment);
                }
            }
        }
    }
    std::cout << "contact segements:" << std::endl;
    for(auto seg: master_segments)
    {
        std::cout << "first:" << seg.first->id << std::endl;
        std::cout << "second:" << seg.second->id << std::endl;
        std::cout << "---" << std::endl;
    }
}

bool  Contact::is_penetrating(const Set<Node*>*, const Node*)
{
    // slave node X is projected onto the master segnment
    // 
    //      __________x0__________ slave segment
    //                |
    //                | e_n is unit vector from master segment intersection to slave segment
    //                |
    //      x1 _______|_________x2  master segment

    // X_c is the point where the node is projected onto the surface.
    // Get search distance from all elements connected to node
    // g = (x - x1)*e_n >= 0
    // Slave node x is projected onto the piecewise linear segments of the
    // master segment with xc as the projected point

    float g = 0;





    // g > 0: no contact
    // g < 0: contact
    return g > 0? false: true;
    
    
}