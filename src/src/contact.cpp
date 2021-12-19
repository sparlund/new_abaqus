#include "../include/contact.h"
#include "../include/element.h"

Contact::Contact(Set<Node*>& master, Set<Node*>& slave) : master(master), slave(slave)
{
    std::cout << "*CONTACT PAIR: slave set = " << slave.get_set_name() << ", master set = " << master.get_set_name() << std::endl;
    // Create segments from nodes in master set
    // not needed for slave set?
    for(size_t i = 0; i < master.size(); i++)
    {
        auto node = master.get_entity(i);
        auto connected_elements = node->connected_elements;
        for(auto element : connected_elements)
        {
            std::cout << element->id << std::endl;
            auto segments = element->get_segments(node);
            for(const auto segment : segments)
            {
                if(master.is_entity_in_set(segment.first) && master.is_entity_in_set(segment.second))
                {
                    // Check we haven't already added this segment, in any direction
                    if(!is_segment_in_master_segments(segment))
                    {
                        master_segments.push_back(segment);
                    }
                }
            }
        }
    }
}

bool Contact::is_segment_in_master_segments(const Segment& segment) const
{
    for(const auto& master_segment: master_segments)
    {
        if((segment.first == master_segment.first && segment.second == master_segment.second) ||
           (segment.second == master_segment.first && segment.first == master_segment.second))
        {
            return true;
        }
    }
    return false;
}

bool  Contact::is_penetrating(const Segment& segment, const Node* node)
{
    // slave node X is projected onto the master segnment
    // don't use concept of segment on slave side
    //
    //      __________x0__________ slave node
    //                |
    //                | e_n is unit vector from master segment intersection to slave segment
    //                |
    //      x1 _______|_________x2  master segment

    // X_c is the point where the node is projected onto the surface.
    // Get search distance from all elements connected to node
    // g = (x2 - x1)*e_n >= 0
    // Slave node x is projected onto the piecewise linear segments of the
    // master segment with xc as the projected point

    // hur vet man om e_n pekar utåt eller inåt?!
    Eigen::Matrix<float,1,2> x1, x2, e_n;
    x1  << segment.first->x,segment.first->y;
    x2  << segment.second->x,segment.second->y;
    e_n << -1*(x2[1] - x2[0]), (x1[1] - x1[0]);
    float g = (x2 - x1)*e_n.transpose();
    auto contact_status = g > 0? false: true;
    std::cout << "Node " << node->id << ": g = " << g << ", contact = " << contact_status << std::endl;
    // g > 0: no contact
    // g < 0: contact
    return contact_status;
    
    
}