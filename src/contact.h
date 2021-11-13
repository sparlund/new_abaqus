#pragma once

#include "node.h"
#include "set.h"

#include <vector>
#include <utility>

// 1) Uses the penalty method
// 2) Slave cannot penetrate master
// 3) Node to segment for 2D contact
// 4) Contact force is an internal force at the interface
using Segment = std::pair<Node&, Node&>;
class Contact
{
private:
    // penalty value is dependent on stiffness of contacting materials! neets to be set dynamically.
    // typically ~1-50*E
    float penalty = 0;
    float allowed penetration = 1e-2;
    size_t number_of_time_steps = 100;
    const Set<Node>& master, slave;
    const Element_dimensions element_dimensions;
    // is the given node penetrating any 2D segment?
    bool  is_penetrating(const Set<Node>&, const Node&);
    // how large is the penetration?
    float gap(const Segment&, const Node&);
    std::vector<Node&> old_nodes_in_contact, current_nodes_in_contact;
    std::vector<Segment> master_segments;
public:
    Contact(Set<Node>& master, Set<Node>& slave);
};


