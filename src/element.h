#pragma once
#include <vector>
#include <memory>
#include "node.h"
#include "pid.h"

// virtual class
class Element{
private:
    
public:
    unsigned int id;
    std::vector<std::shared_ptr<Node>> connectivity;
    std::shared_ptr<Pid> pid;
    virtual std::vector<std::shared_ptr<Node>> get_connectivity()=0;
    virtual std::shared_ptr<Pid> get_pid()=0;
    virtual ~Element(){};
    Element(){};
};
