#pragma once

#include <array>
#include <vector>
#include <algorithm>

#include <line.hpp>
#include <lite_tetrahedron.hpp>

class plane {
public:
    plane() = default;
    plane(size_t res_x, size_t res_y, const std::array<double, 4>& global_boundaries);

    size_t count_all_intersections();

    /*
     * return amount of lines to intersect polygons
     */
    size_t find_intersections_with_tetrahedron(const lite_tetrahedron& tetra, size_t id);

private:
    /*
     * all line equations is on X,Y plane
     */

    static inline double line_common_eq(std::array<double, 3> p1, std::array<double, 3> p2, std::array<double, 3> pos);

    static inline double line_rev_function_eq(std::array<double, 3> p1, std::array<double, 3> p2, double y);

    size_t find_intersections_with_polygon(std::array<std::array<double, 3>, 3> points, size_t id, size_t polygon_id);

    /*
     * [] sequence is [x][y]
     */
    std::vector<std::vector<line>> _lines{};
    std::array<double, 4> _global_boundaries{};

    double _step_x{}, _step_y{};
};