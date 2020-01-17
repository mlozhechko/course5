#pragma once

#include <algorithm>
#include <array>
#include <thread>
#include <vector>
#include <list>

#include <omp.h>

#include <line.hpp>
#include <tetra.hpp>
#include <object3d_base.hpp>
#include <object2d.hpp>

class plane {
public:
    plane() = delete;
    explicit plane(size_t res_x, size_t res_y, std::vector<object3d_base> objects3d, std::vector<double> global_boundaries = {});
    void find_intersections();

    /*
     * trace_all rays and write result to object2d
     */

    object2d trace_rays(tetra_value value_alpha, tetra_value value_Q);

    size_t count_all_intersections();

    size_t get_x() const;
    size_t get_y() const;

    /*
     * debug features:
     */
    void print_all_lines_with_intersection();

private:
    double get_pixel_by_x(double x);
    double get_pixel_by_y(double y);

    /*
     * return amount of lines to intersect polygons
     */
    size_t find_intersections_with_tetrahedron(const tetra& tetra, size_t id,
                                               int internal_thread_id);

    /*
     * all line equations is on X,Y plane
     */
    static double line_common_eq(const double* p1, const double* p2, const double* pos);
    static double line_rev_function_eq(const double* p1, const double* p2, double y);

    size_t find_intersections_with_polygon(std::array<const double *, 3> points, size_t id, size_t polygon_id,
                                           int internal_thread_id, bool is_solid = false, double solid_color = 0);

    std::shared_ptr<std::vector<tetra>> _data{};

    /*
     * lines indexes order is [x][y]
     */
    std::vector<std::vector<line>> _lines{};
    std::array<double, 4> _global_boundaries{};

    size_t _x{}, _y{};

    double _step_x{}, _step_y{};
};