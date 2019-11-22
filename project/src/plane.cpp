#include <plane.hpp>

plane::plane(size_t res_x, size_t res_y, const std::array<double, 4>& global_boundaries)
: _global_boundaries(global_boundaries) {
    double delta_x(global_boundaries[0] - global_boundaries[1]);
    double delta_y(global_boundaries[2] - global_boundaries[3]);

    _step_x = delta_x / (res_x - 1.);
    _step_y = delta_y / (res_y - 1.);

    _lines.resize(res_x);
    double current_x = global_boundaries[1]; //x_min
    for (size_t i = 0; i < res_x; i++) {
        double current_y = global_boundaries[3]; //y_min
        _lines.at(i).reserve(res_y);
        for (size_t j = 0; j < res_y; j++) {
            _lines.at(i).push_back(line(current_x, current_y));
            current_y = current_y + _step_y;
        }
        current_x = current_x + _step_x;
    }
}

size_t plane::count_all_intersections() {
    size_t sum = 0;
    for (auto &i: _lines) {
        for (auto &j: i) {
            sum += j.number_of_intersections();
        }
    }

    return sum;
}

size_t plane::find_intersections_with_tetrahedron(const lite_tetrahedron& tetra, size_t id) {
    /*
     * 0001 0, 1, 2 points
     * 0010 0, 1, 3 points
     * 0100 0, 2, 3 points
     * 1000 1, 2, 3 points
     */

    size_t sum(0);
    sum += find_intersections_with_polygon({tetra[0], tetra[1], tetra[2]}, id, 0);
    sum += find_intersections_with_polygon({tetra[0], tetra[1], tetra[3]}, id, 1);
    sum += find_intersections_with_polygon({tetra[0], tetra[2], tetra[3]}, id, 2);
    sum += find_intersections_with_polygon({tetra[1], tetra[2], tetra[3]}, id, 3);

    if (sum % 2 == 1) {
        throw std::runtime_error("critical error. odd number of intersections");
    }

    return sum / 2;
}

inline double plane::line_common_eq(std::array<double, 3> p1, std::array<double, 3> p2, std::array<double, 3> pos) {
    return (p2[1] - p1[1]) * pos[0] + (p1[0] - p2[0]) * pos[1] + (p2[0] * p1[1] - p1[0] * p2[1]);
}

inline double plane::line_rev_function_eq(std::array<double, 3> p1, std::array<double, 3> p2, double y) {
    return (p1[0] - p2[0]) * (y - p1[1]) / (p1[1] - p2[1]) + p1[0];
}

size_t plane::find_intersections_with_polygon(std::array<std::array<double, 3>, 3> points, size_t id, size_t polygon_id) {
    size_t counter(0);

    std::sort(points.begin(), points.end(), [](auto a, auto b){
        return a[1] > b[1];
    });

    /*
     * find p[1] pos relative to p[0] to p[2] line
     */
    double rel_qual = line_common_eq(points[0], points[2], points[1]);

    /*
     * let L be line passing through points p[0] and p[2]
     *
     * check if L is ascending and p[1] is below the line
     */
    bool asc_above = (points[0][0] >= points[2][0]) && (rel_qual >= 0);

    /*
     * check if L is descending and p[1] is above the line
     */
    bool des_below = (points[0][0] < points[2][0]) && (rel_qual > 0);

    /*
     * p[1] may on the left side or on right side of the line
     * also p[1].y is greater than p[2].y and lower than p[0].y
     *
     * therefore we can fully describe position of triangle of plane, and use it
     * to find grid area that intersects polygon_id.
     *
     * position is true if p[2] is on the right side of L
     */
    bool position = !(asc_above || des_below);

    double y_max = points[0][1];
    double y_min = points[2][1];

    size_t y_max_index = std::floor((y_max - _global_boundaries[3]) / _step_y);
    size_t y_min_index = std::ceil((y_min - _global_boundaries[3]) / _step_y);

    size_t y_index_it = y_min_index;
    double y_it = _lines[0][y_index_it].y();

    /*
     * iterate over triangle's edges and find intersection points
     */
    for (; y_index_it <= y_max_index; y_index_it++) {
        double x_max(0), x_min(0);
        if (position) {
            x_min = line_rev_function_eq(points[0], points[2], y_it);
            if (y_it < points[1][1]) {
                x_max = line_rev_function_eq(points[2], points[1], y_it);
            } else {
                x_max = line_rev_function_eq(points[0], points[1], y_it);
            }
        } else {
            x_max = line_rev_function_eq(points[0], points[2], y_it);
            if (y_it < points[1][1]) {
                x_min = line_rev_function_eq(points[2], points[1], y_it);
            } else {
                x_min = line_rev_function_eq(points[0], points[1], y_it);
            }
        }

        size_t x_max_index = std::floor((x_max - _global_boundaries[1]) / _step_x);
        size_t x_min_index = std::ceil((x_min - _global_boundaries[1]) / _step_x);

        for (size_t i = x_min_index; i <= x_max_index; i++) {
            _lines[i][y_index_it].add_tetra_intersection(id, polygon_id);
            counter++;
        }

        y_it = y_it + _step_y;
    }

    return counter;
}

void plane::trace_rays(const std::vector<lite_tetrahedron>& tetra_vec) {
//    for (size_t i = 0; i < _lines.size(); i++) {
//        for (size_t j = 0; j < _lines[i].size(); j++) {
//            _lines[i][j].calculate_intersections(tetra_vec);
//
//
//            std::cout << "alpha" << _lines[i][j].calculate_ray_value(tetra_vec, tetra_value::alpha) << std::endl;
//        }
//    }

    _lines[2][2].calculate_intersections(tetra_vec);
    std::cout << "alpha" << _lines[2][2].calculate_ray_value(tetra_vec, tetra_value::alpha) << std::endl;
}
