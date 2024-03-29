#include "../include/contact.h"
#include "../include/element.h"

#include <iomanip>


inline bool SameSign(double a, double b) {
    return a*b >= 0.d;
}


std::vector<std::pair<Node*, double>> Contact::get_penetrating_nodes()
{
    // Check if any slave node is penetrating any master set segments
    // and correct position of those who do
    for (auto& node: slave->entities)
    {
        for (const auto& segment : master_segments)
        {
            auto g = gap(segment,node);
            // g > 0: no contact
            // g < 0: contact
            // Check if node is on the same side of the segment as initally.
            // Has it changed side it's penetrating
            auto current_distance_sign  = ((segment.second->x - segment.first->x)*(node->y - segment.first->y) - (segment.second->y - segment.first->y)*(node->x - segment.first->x));
            auto original_distance_sign = ((segment.second->x - segment.first->x)*(node->original_y - segment.first->y) - (segment.second->y - segment.first->y)*(node->original_x - segment.first->x));
            auto penetration = SameSign(current_distance_sign,original_distance_sign)? false: true;
            if (penetration)
            {

            }
        }
    }
    
}


Contact::Contact(Set<Node*>* master_, Set<Node*>* slave_) : master(master_), slave(slave_)
{
    std::cout << "*CONTACT PAIR: slave set = " << slave->name << ", master set = " << master->name << std::endl;
    // Create segments from nodes in master set
    // not needed for slave set, there we go node by node
    auto node = slave->entities[0];
    for(size_t i = 0; i < master->size(); i++)
    {
        auto* node = master->get_entity(i);
        auto connected_elements = node->connected_elements;
        for(const auto& element : connected_elements)
        {
            auto segments = element->get_segments(node);
            for(const auto segment : segments)
            {
                if(master->is_entity_in_set(segment.first) && master->is_entity_in_set(segment.second))
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

double Contact::gap(const Segment& segment, const Node* node)
{
    // https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
    // https://stackoverflow.com/questions/849211/short
    double x = node->x;
    double y = node->y;
    double x1 = segment.first->x;
    double y1 = segment.first->y;
    double x2 = segment.second->x;
    double y2 = segment.second->y;
    double A = x - x1;
    double B = y - y1;
    double C = x2 - x1;
    double D = y2 - y1;
    double dot = (A * C) + (B * D);
    double len_sq = (C * C) + (D * D);
    double param = -1;
    // in case of 0 length line
    if (len_sq != 0)
    {
        param = dot / len_sq;
    }
    double xx, yy;

    if (param < 0)
    {
        xx = x1;
        yy = y1;
    }
    else if (param > 1)
    {
        xx = x2;
        yy = y2;
    }
    else
    {
        xx = x1 + (param * C);
        yy = y1 + (param * D);
    }
    double dx = x - xx;
    double dy = y - yy;
    return std::sqrt((dx * dx) + (dy * dy));
}