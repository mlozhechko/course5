#pragma once

#include <array>
#include <vector>
#include <algorithm>
#include <thread>

#include <omp.h>

#include <line.hpp>
#include <tetra.hpp>

using float_matrix = std::vector<std::vector<float>>;

class plane {
public:
    plane() = default;
    plane(size_t res_x, size_t res_y, const std::array<double, 4>& global_boundaries);

    size_t count_all_intersections();

    void find_intersections_with_tetra_vector(const std::vector<tetra>& tetra_vec);

    /*
     * trace_all rays
     *
     * TODO:
     * 1. incapsulate tetra vector
     */
    std::pair<float_matrix, float_matrix>
    trace_rays(const std::vector<tetra>& tetra_vec, tetra_value value_alpha,
               tetra_value value_Q);

    size_t get_x() {
        return _x;
    }

    size_t get_y() {
        return _y;
    }

    /*
     * debug
     */
    void print_all_lines_with_intersection();
private:
    /*
     * return amount of lines to intersect polygons
     */
    size_t find_intersections_with_tetrahedron(const tetra& tetra, size_t id, int internal_thread_id);

    /*
     * all line equations is on X,Y plane
     */
    static double line_common_eq(const double *p1, const double *p2, const double *pos);
    static double line_rev_function_eq(const double *p1, const double *p2, double y);

    size_t find_intersections_with_polygon(std::array<const double *, 3> points, size_t id, size_t polygon_id,
                                           int internal_thread_id);

    /*
     * [] sequence is [x][y]
     */
    std::vector<std::vector<line>> _lines{};
    std::array<double, 4> _global_boundaries{};

    size_t _x{}, _y{};

    double _step_x{}, _step_y{};
};