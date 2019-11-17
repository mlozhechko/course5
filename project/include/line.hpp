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

    void calculate_intersections_alpha()

private:
    double _x, _y;
    bool alternation_flag{false};

    /*
     * std::bitset<32> structure:
     * 4 bits | 28 bits
     * first 4 bits is polygon intersection info
     * there are 4 polygons in tetrahedron, n bit is set if line intersects n-th polygon if tetrahedron
     *
     * 28 last bits is id number of tetrahedron
     */
    std::vector<std::bitset<32>> _tetra_intersections;
};