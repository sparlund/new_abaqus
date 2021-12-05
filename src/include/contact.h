#pragma once

#include "node.h"
#include "set.h"

#include <vector>
#include <utility>

// 1) Uses the penalty method
// 2) Slave cannot penetrate master
// 3) Node to segment for 2D contact
// 4) Contact force is an internal force at the interface
using Segment = std::pair<Node*, Node*>;
class Contact
{
private:
    // penalty value is dependent on stiffness of contacting materials! neets to be set dynamically.
    // typically ~1-50*E
    float penalty = 0;
    // tolerance for equilibrium iterations. Maybe update dynamically?
    float relative_tolerance = 1e-6;
    // max allowed penetration. too large obviously bad, but too small and difficult its to solve
    float allowed_penetration = 1e-2;
    size_t number_of_time_steps = 100;
    // is the given node penetrating any 2D segment?
    bool  is_penetrating(const Set<Node*>*, const Node*);
    // how large is the penetration?
    float gap(const Segment&, const Node*);
    std::vector<Node*> old_nodes_in_contact, current_nodes_in_contact;
    std::vector<Segment> master_segments;
public:
    const Set<Node*>& master, slave;
    Contact(Set<Node*>& master, Set<Node*>& slave);
};


