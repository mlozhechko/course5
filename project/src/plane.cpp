
#include <plane.hpp>
#include <atomic>

plane::plane(size_t res_x, size_t res_y, const std::array<double, 4>& global_boundaries)
    : _global_boundaries(global_boundaries) {
    _x = res_x;
    _y = res_y;

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
    for (auto& i: _lines) {
        for (auto& j: i) {
            sum += j.number_of_intersections();
        }
    }

    return sum;
}

size_t plane::find_intersections_with_tetrahedron(const tetra& tetra, size_t id, int internal_thread_id) {
    /*
     * 0001 0, 1, 2 points
     * 0010 0, 1, 3 points
     * 0100 0, 2, 3 points
     * 1000 1, 2, 3 points
     */

    size_t sum(0);
    sum += find_intersections_with_polygon({tetra[0].data(), tetra[1].data(), tetra[2].data()}, id, 0, internal_thread_id);
    sum += find_intersections_with_polygon({tetra[0].data(), tetra[1].data(), tetra[3].data()}, id, 1, internal_thread_id);
    sum += find_intersections_with_polygon({tetra[0].data(), tetra[2].data(), tetra[3].data()}, id, 2, internal_thread_id);
    sum += find_intersections_with_polygon({tetra[1].data(), tetra[2].data(), tetra[3].data()}, id, 3, internal_thread_id);

    if (sum % 2 == 1) {
        throw std::runtime_error("critical error. odd number of intersections");
    }

    return sum / 2;
}

double plane::line_common_eq(const double *p1, const double *p2, const double *pos) {
    return (p2[1] - p1[1]) * pos[0] + (p1[0] - p2[0]) * pos[1] + (p2[0] * p1[1] - p1[0] * p2[1]);
}

double plane::line_rev_function_eq(const double *p1, const double *p2, double y) {
    return (p1[0] - p2[0]) * (y - p1[1]) / (p1[1] - p2[1]) + p1[0];
}

size_t plane::find_intersections_with_polygon(std::array<const double *, 3> points, size_t id, size_t polygon_id,
                                              int internal_thread_id) {
    size_t counter(0);

    std::sort(points.begin(), points.end(), [](auto& a, auto& b) {
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

    const size_t y_max_index = std::floor((y_max - _global_boundaries[3]) / _step_y);
    const size_t y_min_index = std::ceil((y_min - _global_boundaries[3]) / _step_y);

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

        const size_t x_max_index = std::floor((x_max - _global_boundaries[1]) / _step_x);
        const size_t x_min_index = std::ceil((x_min - _global_boundaries[1]) / _step_x);

        for (size_t i = x_min_index; i <= x_max_index; i++) {
            _lines[i][y_index_it].add_tetra_intersection(id, polygon_id, internal_thread_id);
            counter++;
        }

        y_it = y_it + _step_y;
    }

    return counter;
}

std::pair<float_matrix, float_matrix>
plane::trace_rays(const std::vector<tetra>& tetra_vec, const tetra_value value_alpha,
                  const tetra_value value_Q) {
    std::vector<std::vector<float>> result_x{};
    std::vector<std::vector<float>> result_y{};
    result_x.resize(_lines.size());
    result_y.resize(_lines.size());

    const size_t lines_size_1d = _lines.size();
    if (lines_size_1d == 0) {
        throw std::runtime_error("critical error. empty plane");
    }
    const size_t lines_size_2d = _lines[0].size();

    for (size_t i = 0; i < lines_size_1d; i++) {
        result_x[i].resize(lines_size_2d);
        result_y[i].resize(lines_size_2d);
    }

#pragma omp parallel for default(none) shared(result_x, result_y, tetra_vec) num_threads(AMOUNT_OF_THREADS) schedule(dynamic, 8) collapse(2)
    for (size_t i = 0; i < lines_size_1d; i++) {
        for (size_t j = 0; j < lines_size_2d; j++) {
            _lines[i][j].calculate_intersections(tetra_vec);
            result_x[i][j] = _lines[i][j].direct_calculate_ray_value(tetra_vec, value_alpha);
            result_y[i][j] = _lines[i][j].integrate_ray_value_by_i(tetra_vec, value_alpha, value_Q);

        }
    }

    return {result_x, result_y};
}

void plane::print_all_lines_with_intersection() {
    for (auto& i: _lines) {
        for (auto& j: i) {
            if (j.number_of_intersections() > 0) {
                std::cout << "" << j.x() << " " << j.y() << "" << std::endl;
            }
        }
    }
}

void plane::find_intersections_with_tetra_vector(const std::vector<tetra>& tetra_vec) {
    size_t tetrahedron_vector_size = tetra_vec.size();

#pragma omp parallel for default(none) shared(tetrahedron_vector_size, tetra_vec) num_threads(AMOUNT_OF_THREADS) schedule(dynamic, 8)
    for (size_t i = 0; i < tetrahedron_vector_size; i++) {
        int tid = omp_get_thread_num();
        find_intersections_with_tetrahedron(tetra_vec[i], i, tid);
    }
}