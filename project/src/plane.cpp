
#include <atomic>
#include <plane.hpp>

//plane::plane(size_t res_x, size_t res_y, const std::array<double, 4>& global_boundaries)
//    : _global_boundaries(global_boundaries) {
//    _x = res_x;
//    _y = res_y;
//
//    double delta_x(global_boundaries[0] - global_boundaries[1]);
//    double delta_y(global_boundaries[2] - global_boundaries[3]);
//
//    _step_x = delta_x / (res_x - 1.);
//    _step_y = delta_y / (res_y - 1.);
//
//    _lines.resize(res_x);
//    double current_x = global_boundaries[1]; // x_min
//    for (size_t i = 0; i < res_x; i++) {
//        double current_y = global_boundaries[3]; // y_min
//        _lines.at(i).reserve(res_y);
//        for (size_t j = 0; j < res_y; j++) {
//            _lines.at(i).push_back(line(current_x, current_y));
//            current_y = current_y + _step_y;
//        }
//        current_x = current_x + _step_x;
//    }
//}

size_t plane::count_all_intersections() {
    size_t sum = 0;
    for (auto& i : _lines) {
        for (auto& j : i) {
            sum += j.number_of_intersections();
        }
    }

    return sum;
}

size_t plane::find_intersections_with_tetrahedron(const tetra& tetra, size_t id,
                                                  int internal_thread_id) {
    /*
     * 0001 0, 1, 2 points
     * 0010 0, 1, 3 points
     * 0100 0, 2, 3 points
     * 1000 1, 2, 3 points
     */

    size_t sum(0);
    sum += find_intersections_with_polygon({tetra[0].data(), tetra[1].data(), tetra[2].data()}, id,
                                           0, internal_thread_id);
    sum += find_intersections_with_polygon({tetra[0].data(), tetra[1].data(), tetra[3].data()}, id,
                                           1, internal_thread_id);
    sum += find_intersections_with_polygon({tetra[0].data(), tetra[2].data(), tetra[3].data()}, id,
                                           2, internal_thread_id);
    sum += find_intersections_with_polygon({tetra[1].data(), tetra[2].data(), tetra[3].data()}, id,
                                           3, internal_thread_id);

    if (sum % 2 == 1) {
        throw std::runtime_error("critical error. odd number of intersections");
    }

    return sum / 2;
}

double plane::line_common_eq(const double* p1, const double* p2, const double* pos) {
    return (p2[1] - p1[1]) * pos[0] + (p1[0] - p2[0]) * pos[1] + (p2[0] * p1[1] - p1[0] * p2[1]);
}

double plane::line_rev_function_eq(const double* p1, const double* p2, double y) {
    return (p1[0] - p2[0]) * (y - p1[1]) / (p1[1] - p2[1]) + p1[0];
}

size_t plane::find_intersections_with_polygon(std::array<const double*, 3> points, size_t id,
                                              size_t polygon_id, int internal_thread_id) {
    size_t counter(0);

    std::sort(points.begin(), points.end(), [](auto& a, auto& b) { return a[1] > b[1]; });

    /*
     * find p[1] pos relative to p[0] to p[2] line
     */
    double rel_qual = line_common_eq(points[0], points[2], points[1]);

    /*
     * let L be line passing through points p[0] and p[2]
     *
     * check if L is ascending and p[1] is below the line
     */
    bool asc_above = ((points[0][0] >= points[2][0]) && (rel_qual >= 0));

    /*
     * check if L is descending and p[1] is above the line
     */
    bool des_below = ((points[0][0] < points[2][0]) && (rel_qual > 0));

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

    //    const size_t y_max_index = std::floor((y_max - _global_boundaries[3]) / _step_y);
    //    const size_t y_min_index = std::ceil((y_min - _global_boundaries[3]) / _step_y);
    const size_t y_max_index = std::floor(get_pixel_by_y(y_max));
    const size_t y_min_index = std::ceil(get_pixel_by_y(y_min));

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
        //
        //        const size_t x_max_index = std::floor((x_max - _global_boundaries[1]) / _step_x);
        //        const size_t x_min_index = std::ceil((x_min - _global_boundaries[1]) / _step_x);
        const size_t x_max_index = std::floor(get_pixel_by_x(x_max));
        const size_t x_min_index = std::ceil(get_pixel_by_x(x_min));

        for (size_t i = x_min_index; i <= x_max_index; i++) {
            _lines[i][y_index_it].add_tetra_intersection(id, polygon_id, internal_thread_id);
            counter++;
        }

        y_it = y_it + _step_y;
    }

    return counter;
}

object2d plane::trace_rays(const tetra_value value_alpha, const tetra_value value_Q) {
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

#pragma omp parallel for default(none) shared(result_x, result_y, _data)                       \
    num_threads(AMOUNT_OF_THREADS) schedule(dynamic, 8) collapse(2)
    for (size_t i = 0; i < lines_size_1d; i++) {
        for (size_t j = 0; j < lines_size_2d; j++) {
            _lines[i][j].calculate_intersections(*_data);
            result_x[i][j] = _lines[i][j].direct_calculate_ray_value(*_data, value_alpha);
            result_y[i][j] = _lines[i][j].integrate_ray_value_by_i(*_data, value_alpha, value_Q);
            _lines[i][j].free_memory();
        }
    }

    return object2d{std::pair{result_x, result_y}};
}

void plane::print_all_lines_with_intersection() {
    for (auto& i : _lines) {
        for (auto& j : i) {
            if (j.number_of_intersections() > 0) {
                std::cout << "" << j.x() << " " << j.y() << "" << std::endl;
            }
        }
    }
}

void plane::find_intersections() {
    size_t tetrahedron_vector_size = _data->size();

#pragma omp parallel for default(none) shared(tetrahedron_vector_size, _data)                  \
    num_threads(AMOUNT_OF_THREADS) schedule(dynamic, 8)
    for (size_t i = 0; i < tetrahedron_vector_size; i++) {
        int tid = omp_get_thread_num();
        find_intersections_with_tetrahedron((*_data)[i], i, tid);
    }
}

/*
 * accretion disk is sphere
 *
 * equation: (x - x0)^2 + (y - y0)^2 == acc_r^2
 */
//void plane::mark_accretion_disk_space() {
//    const size_t y_max_index = std::floor(get_pixel_by_y(acc_y0 + acc_r));
//    const size_t y_min_index = std::ceil(get_pixel_by_y(acc_y0 - acc_r));
//    double current_y_it = acc_y0 - acc_r;
//
//    for (size_t y_it = y_min_index; y_it <= y_max_index; y_it++) {
//        /* xp = (x - x0) = sqrt(acc_r^2 - (y - y0)^2) */
//        double yp = (current_y_it - acc_y0);
//        double xp = sqrt(acc_r * acc_r - yp * yp);
//
//        double x_max_it = xp + acc_x0;
//        double x_min_it = acc_x0 - xp;
//
//        const size_t x_max_index = std::floor(get_pixel_by_x(x_max_it));
//        const size_t x_min_index = std::ceil(get_pixel_by_x(x_min_it));
//
//        for (size_t x_it = x_min_index; x_it <= x_max_index; x_it++) {
//            _lines[x_it][y_it].mark_void(marked_accretor_value);
//        }
//
//        current_y_it += _step_y;
//    }
//}

double plane::get_pixel_by_x(double x) {
    return (x - _global_boundaries[1]) / _step_x;
}

double plane::get_pixel_by_y(double y) {
    return (y - _global_boundaries[3]) / _step_y;
}

double vector_2_norm(const std::array<double, 3>& vec) {
    double sum{};
    for (const auto& i : vec) {
        sum += i * i;
    }
    return sqrt(sum);
}

std::array<double, 3> vector_multiplication(const std::array<double, 3>& m1,
                                            const std::array<double, 3>& m2) {
    std::array<double, 3> res{};
    res[0] = m1[1] * m2[2] - m1[2] * m2[1];
    res[1] = -(m1[0] * m2[2]) + (m1[2] * m2[0]);
    res[2] = m1[0] * m2[1] - m1[1] * m2[0];

    return res;
}

double roche_lobe_potential(const std::array<double, 3>& r, double acc_x, double donor_x,
                            double mc_x) {
    double acc_denominator = vector_2_norm({r[0] - acc_x, r[1], r[2]});
    //    std::cout << "acc denominator: " << acc_denominator << std::endl;
    double donor_denominator = vector_2_norm({r[0] - donor_x, r[1], r[2]});
    //    std::cout << "donor denominator: " << donor_denominator << std::endl;
    double tmp_omega_mult =
        vector_2_norm(vector_multiplication({r[0] - mc_x, r[1], r[2]}, {0, OMEGA, 0}));
    double omega = (1. / 2.) * tmp_omega_mult * tmp_omega_mult;
    //    std::cout << "omega: " << omega << std::endl;

    double F = -((G_SOL * Ma) / acc_denominator) - ((G_SOL * Md) / donor_denominator) - omega;
    return F;
}

std::array<double, 3> rotate_vector_by_y_axis(const std::array<double, 3>& line_vector,
                                              double alpha) {
    std::array<double, 3> res{};
    //    double tmp = line_vector[0];
    res[0] = line_vector[0] * cos(alpha) + line_vector[2] * sin(alpha);
    res[1] = line_vector[1];
    res[2] = -line_vector[0] * sin(alpha) + line_vector[2] * cos(alpha);

    return res;
}

std::array<double, 3> rotate_vector_by_z_axis(const std::array<double, 3>& line_vector,
                                              double alpha) {
    std::array<double, 3> res{};
    //    double tmp = line_vector[0];
    res[0] = line_vector[0] * cos(alpha) + line_vector[1] * sin(alpha);
    res[1] = -line_vector[0] * sin(alpha) + line_vector[1] * cos(alpha);
    res[2] = line_vector[2];

    return res;
}

void add_vector(std::array<double, 3>& arr1, const std::array<double, 3>& arr2) {
    for (size_t i = 0; i < arr1.size(); i++) {
        arr1[i] += arr2[i];
    }
}

void print_vector(const std::array<double, 3>& vec) {
    std::cout << "(" << vec[0] << "," /* << vec[1] << " "*/ << vec[2] << ")" << std::endl;
}
//
//
//std::vector<tetra> plane::build3d_model_donor_roche_lobe() {
//    /*
//     * assume that in the beginning roche lobe placed on x axis
//     *
//     * build 3d model and then rotate it
//     */
//
//    double donor_pos_x = 1 - L;
//    double mass_center_pos_x = (donor_pos_x * Md + acc_x0 * Ma) / (Ma + Md);
//    double alpha = Md / (Ma + Md);
//    //    double lagrange1_pos_x = mass_center_pos_x - L * (1 - pow((alpha / 3.), 1. / 3.));
//    //    std::cout << lagrange1_pos_x << std::endl;
//
//    /*
//     * for some reason defined lagrange_pos_x is another value
//     */
//    double lagrange1_pos_x = 0.35515;
//
//    /*
//     * let norm multiply value by G_SOL
//     * by definition OMEGA is perpendicular to r
//     *
//     * using formula (4) https://arxiv.org/pdf/1702.00587.pdf
//     */
//
//    /*
//     * let's find potential
//     */
//    const double l1_potential =
//        roche_lobe_potential({lagrange1_pos_x, 0, 0}, acc_x0, donor_pos_x, mass_center_pos_x);
//    //    std::cout << l1_potential << std::endl;
//
//    const std::array<double, 3> def_step_vector{0.001, 0, 0};
//
//    double step_angle_x = PI / 256;
//    double step_angle_y = PI / 256;
//
//    double angle_x = 0;
//    double angle_y = -PI + step_angle_y;
//
//    std::array<double, 3> bottom_point{};
//    std::array<double, 3> top_point{};
//
//    {
//        std::array<double, 3> point_trace_vec{donor_pos_x, 0, 0};
//        std::array<double, 3> step_vector{0, 0, 0.001};
//
//        double f_value{};
//        do {
//            add_vector(point_trace_vec, step_vector);
//            f_value = roche_lobe_potential(point_trace_vec, acc_x0, donor_pos_x, mass_center_pos_x);
//        } while (f_value < l1_potential);
//
//        top_point = point_trace_vec;
//        point_trace_vec = {donor_pos_x, 0, 0};
//        step_vector = {0, 0, -0.001};
//
//        do {
//            add_vector(point_trace_vec, step_vector);
//            f_value = roche_lobe_potential(point_trace_vec, acc_x0, donor_pos_x, mass_center_pos_x);
//        } while (f_value < l1_potential);
//
//        bottom_point = point_trace_vec;
//    }
//
//    std::vector<std::vector<std::array<double, 3>>> lobe_points{};
//
//    /*
//     * TODO:
//     * 1. resize vectors in tracing
//     * 2. resize lobe 3d model
//     * 3. parallelize algorithm
//     * 4. refactor code (!)
//     * 5. refactor tetra scalar values, implement base class that don't explicitly rely on q and
//     * alpha, implement derivative classes for roche_lobe and star center
//     * 6. implement write 3d model object data back to .vtk (unstructured_grid_tetra)
//     */
//
//    while (angle_y < (PI - step_angle_y + std::numeric_limits<double>::epsilon())) {
//        std::array<double, 3> y_step_vector = rotate_vector_by_z_axis(def_step_vector, angle_y);
//        angle_y += step_angle_y;
//
//        std::vector<std::array<double, 3>> line{};
//        while (angle_x < 2 * PI - step_angle_x + std::numeric_limits<double>::epsilon()) {
//            std::array<double, 3> point_trace_vec{donor_pos_x, 0, 0};
//            std::array<double, 3> step_vector = rotate_vector_by_y_axis(y_step_vector, angle_x);
//
//            double f_value{};
//            do {
//                add_vector(point_trace_vec, step_vector);
//                f_value =
//                    roche_lobe_potential(point_trace_vec, acc_x0, donor_pos_x, mass_center_pos_x);
//            } while (f_value < l1_potential);
//
//            line.push_back(point_trace_vec);
//            angle_x += step_angle_x;
//        }
//
//        lobe_points.push_back(std::move(line));
//        angle_x = 0;
//    }
//
//    /*
//     * create vector of tetras of lobe
//     *
//     * main cycle from bottom to top
//     */
//
//    const std::array<double, 3> center{donor_pos_x, 0, 0};
//    std::vector<tetra> lobe_3d_model{};
//
//    const size_t line_len = lobe_points[0].size();
//    for (size_t i = 1; i < line_len; i++) {
//        lobe_3d_model.emplace_back(std::array<std::array<double, 3>, 4>{center, bottom_point,
//                                                                        lobe_points[0][i],
//                                                                        lobe_points[0][i - 1]},
//                                   10., 10.);
//    }
//    lobe_3d_model.emplace_back(std::array<std::array<double, 3>, 4>{center, bottom_point,
//                                                                    lobe_points[0][0],
//                                                                    lobe_points[0][line_len - 1]},
//                               10., 10.);
//
//    const size_t h_len = lobe_points.size() - 1;
//    for (size_t i = 1; i < line_len; i++) {
//        lobe_3d_model.emplace_back(std::array<std::array<double, 3>, 4>{center, top_point,
//                                                                        lobe_points[h_len][i],
//                                                                        lobe_points[h_len][i - 1]},
//                                   10., 10.);
//    }
//    lobe_3d_model.emplace_back(std::array<std::array<double, 3>, 4>{center, top_point,
//                                                                    lobe_points[h_len][0],
//                                                                    lobe_points[0][line_len - 1]},
//                               10., 10.);
//
//    for (size_t i = 1; i < (h_len + 1); i++) {
//        for (size_t j = 1; j < line_len; j++) {
//            lobe_3d_model.emplace_back(
//                std::array<std::array<double, 3>, 4>{center, lobe_points[i - 1][j - 1],
//                                                     lobe_points[i - 1][j], lobe_points[i][j - 1]},
//                10., 10.);
//            lobe_3d_model.emplace_back(
//                std::array<std::array<double, 3>, 4>{center, lobe_points[i][j - 1],
//                                                     lobe_points[i][j], lobe_points[i - 1][j]},
//                10., 10.);
//        }
//
//        lobe_3d_model.emplace_back(
//            std::array<std::array<double, 3>, 4>{center, lobe_points[i - 1][line_len - 1],
//                                                 lobe_points[i - 1][0],
//                                                 lobe_points[i][line_len - 1]},
//            10., 10.);
//        lobe_3d_model.emplace_back(
//            std::array<std::array<double, 3>, 4>{center, lobe_points[i][line_len - 1],
//                                                 lobe_points[i][0], lobe_points[i - 1][0]},
//            10., 10.);
//    }
//
//    return lobe_3d_model;
//}


plane::plane(size_t res_x, size_t res_y, std::vector<object3d_base> objects3d) {
    if (objects3d.empty()) {
        _global_boundaries = {1, 0, 1, 0};
        std::cerr << "warning! empty set of objects to render" << std::endl;
    } else {
        _global_boundaries = objects3d.at(0).get_boundaries();
        _data = objects3d.at(0).get_pointer();

        const size_t objects3d_size =  objects3d.size();
        for (size_t i = 1; i < objects3d_size; i++) {
            auto tmp = objects3d.at(i).get_boundaries();
            _global_boundaries[0] = std::max(_global_boundaries[0], tmp[0]);
            _global_boundaries[1] = std::min(_global_boundaries[1], tmp[1]);
            _global_boundaries[2] = std::max(_global_boundaries[2], tmp[2]);
            _global_boundaries[3] = std::min(_global_boundaries[3], tmp[3]);

            auto tmp_data = objects3d.at(i).get_pointer();
            _data->insert(_data->end(), std::make_move_iterator(tmp_data->begin()), std::make_move_iterator(tmp_data->end()));
        }
    }

    _x = res_x;
    _y = res_y;

    double delta_x(_global_boundaries[0] - _global_boundaries[1]);
    double delta_y(_global_boundaries[2] - _global_boundaries[3]);

    _step_x = delta_x / (res_x - 1.);
    _step_y = delta_y / (res_y - 1.);

    _lines.resize(res_x);
    double current_x = _global_boundaries[1]; // x_min
    for (size_t i = 0; i < res_x; i++) {
        double current_y = _global_boundaries[3]; // y_min
        _lines.at(i).reserve(res_y);
        for (size_t j = 0; j < res_y; j++) {
            _lines.at(i).push_back(line(current_x, current_y));
            current_y = current_y + _step_y;
        }
        current_x = current_x + _step_x;
    }
}
