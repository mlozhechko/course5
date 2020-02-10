#pragma once

#include <algorithm>
#include <bitset>
#include <iostream>
#include <limits>
#include <mutex>
#include <vector>

#include <config.hpp>
#include <tetra.hpp>

/*
 * temporary solution
 *
 * add tetra intersection custom vector thread safe container
 */

struct intersection_data {
    double delta_z;
    unsigned tetra_id;
};

class line {
public:
    explicit line(double x, double y);

    line(line&) = delete;
    line(line&&) noexcept;

    double x() const;
    double y() const;

    void mark_solid_color(double mark_value);

    void add_tetra_intersection(size_t id, size_t polygon_id, int internal_thread_id);
    size_t number_of_intersections();

    void calculate_intersections(const std::vector<tetra>&);

    /*
     * TODO:
     * 1. make std::function acceptor to all calculators and integrators
     */

    /*
     * sum delta_z within tetrahedron * alpha for all tetrahedrons which
     * intersects the line
     */
    double direct_calculate_ray_value(const std::vector<tetra>& tetra_vector,
                                      tetra_value value_signature);

    /*
     * solve I' + alpha * I = Q equation
     */
    double integrate_ray_value_by_i(const std::vector<tetra>& tetra_vector,
                                    tetra_value alpha_signature, tetra_value q_signature);

    void free_memory();

private:
    double find_polygon_intersection_z(const std::array<double, 3>& p1,
                                       const std::array<double, 3>& p2,
                                       const std::array<double, 3>& p3);

    double _x, _y;
    void ts_tetra_intersections_pushback(std::bitset<32> data);

    /*
     * std::bitset<32> structure:
     * 4 bits | 28 bits
     * first 4 bits is polygon intersection info
     * there are 4 polygons in tetrahedron, n-th bit is set if line intersects n-th
     * polygon if tetrahedron
     *
     * 28 last bits is id number of tetrahedron
     */

    std::vector<std::bitset<32>> _tetra_intersections;
    std::vector<intersection_data> _intersections_delta;

    std::array<std::bitset<32>, MAX_NUMBER_OF_THREADS> buffer_data{};
    std::array<bool, MAX_NUMBER_OF_THREADS> buffer_flags{};

    std::mutex _tetra_intersections_mutex;

    bool _marked_solid{false};
    double _mark_value{};
};