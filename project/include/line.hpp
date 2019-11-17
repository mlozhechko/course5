#pragma once

#include <vector>
#include <bitset>
#include <iostream>

class line {
public:
    explicit line(double x, double y);

    double x();
    double y();

    void add_tetra_intersection(size_t id, size_t polygon_id);
    size_t number_of_intersections();


private:
    double _x, _y;
    bool alternation_flag{false};

    /*
     *
     */
    std::vector<std::bitset<32>> _tetra_intersections;
};