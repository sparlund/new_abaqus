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
    // typically ~1e-50*E
    double penalty = 0;
    // tolerance for equilibrium iterations. Maybe update dynamically?
    double relative_tolerance = 1e-6;
    // max allowed penetration. too large obviously bad, but too small and difficult its to solve
    double allowed_penetration = 1e-2;
    size_t number_of_time_steps = 100;
    // how large is the penetration?
    double gap(const Segment&, const Node*);
    std::vector<Node*> old_nodes_in_contact, current_nodes_in_contact;
    bool is_segment_in_master_segments(const Segment& segment) const;
public:
    double get_penalty_factor() const {return penalty;}
    std::vector<Segment> master_segments;
    // TODO: this should be something like 1e-5, will change when contact works.
    const double residual_tolerance = 0.1;
    Set<Node*>* master;
    Set<Node*>* slave;
    Contact(Set<Node*>* master, Set<Node*>* slave);
    std::vector<std::pair<Node*, double>> get_penetrating_nodes();
};


