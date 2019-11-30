#pragma once

#include <vector>
#include <bitset>
#include <iostream>
#include <algorithm>
#include <limits>

#include <lite_tetrahedron.hpp>

struct intersection_data {
    double delta_z;
    unsigned tetra_id;
};

class line {
public:
    explicit line(double x, double y);

    double x();
    double y();

    void add_tetra_intersection(size_t id, size_t polygon_id);
    size_t number_of_intersections();

    void calculate_intersections(const std::vector<lite_tetrahedron>&);

    /*
     * sum delta_z within tetrahedron * alpha for all tetrahedrons which intersects the line
     */
    double direct_calculate_ray_value(const std::vector<lite_tetrahedron>& tetra_vector, tetra_value value_signature);

    /*
     * solve I' + alpha * I = Q equation
     */
    double integrate_ray_value_by_i(const std::vector<lite_tetrahedron>& tetra_vector,
                                    tetra_value alpha_signature,
                                    tetra_value q_signature);

private:
    double find_polygon_intersection_z(
        const std::array<double, 3>& p1,
        const std::array<double, 3>& p2,
        const std::array<double, 3>& p3);

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

    /*
     * TODO:
     * 1. concat to one vector
     */
    std::vector<std::bitset<32>> _tetra_intersections;
    std::vector<intersection_data> _intersections_delta;
};