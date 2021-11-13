#pragma once

class Output
{
private:
public:
    enum class nodal
    {
        U, // displacement
    }
    enum class element
    {
        S, // all stess components
    }
    Output(/* args */);
    ~Output();
};

Output::Output(/* args */)
{
}

Output::~Output()
{
}


