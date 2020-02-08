#include <line.hpp>

line::line(double x, double y) : _x(x), _y(y) {
    _tetra_intersections.reserve(50);
    //
    //    for (size_t i = 0; i < AMOUNT_OF_THREADS; i++) {
    //        _threads_buffers[i].reserve(10);
    //    }
}

double line::x() const {
    return _x;
}

double line::y() const {
    return _y;
}

const int DELIMITER_POS = 28;

/*
 * 0001 0, 1, 2 points
 * 0010 0, 1, 3 points
 * 0100 0, 2, 3 points
 * 1000 1, 2, 3 points
 */
static std::bitset<32> mask{0x0FFFFFFF};

void line::add_tetra_intersection(size_t id, size_t polygon_id, int internal_thread_id) {
    if (_marked_solid) {
        return;
    }

    buffer_data[internal_thread_id].set(DELIMITER_POS + polygon_id);
    buffer_data[internal_thread_id] |= id;

    /*
     * check for data races and other bad things
     */
    if (!buffer_flags[internal_thread_id]) {
        if ((buffer_data[internal_thread_id] & mask).to_ulong() != id) {
            std::cout << internal_thread_id << std::endl;
            std::cout << (buffer_data[internal_thread_id] & mask).to_ulong() << " " << id
                      << std::endl;
            throw std::runtime_error("tetrahedron intersection fatal data error");
        }
    }

    if (buffer_flags[internal_thread_id]) {
        ts_tetra_intersections_pushback(buffer_data[internal_thread_id]);
        buffer_data[internal_thread_id] = 0;
    }

    buffer_flags[internal_thread_id] = !buffer_flags[internal_thread_id];

    /*
    if (!alternation_flag) {
        _tetra_intersections.emplace_back(id);
    }

    _tetra_intersections.back().set(28 + polygon_id);
    alternation_flag = !alternation_flag;


    //_cont_tetra_intersections.add_intersection(0, id | (polygon_id << 28u));
    */
}

size_t line::number_of_intersections() {
    return _tetra_intersections.size();
}

struct set_intersection_data {
    double key_z;
    intersection_data data;
};

struct set_intersection_data_cmp {
    bool operator()(const set_intersection_data& a, const set_intersection_data& b) const {
        return a.key_z > b.key_z;
    }
};

void line::calculate_intersections(const std::vector<tetra>& tetra_vector) {
    //    boost::container::flat_set<set_intersection_data,
    //    set_intersection_data_cmp> values; std::cout << "!" <<
    //    tetra_vector.size() << std::endl;

    std::vector<set_intersection_data> values;

    size_t intrs_size = _tetra_intersections.size();
    values.reserve(intrs_size);

    for (const auto& it : _tetra_intersections) {
        size_t tetra_id = (it & mask).to_ulong();

        //        std::cout << tetra_id << std::endl;

        size_t i = 0;
        std::array<double, 2> points{};

        const auto& tetra = tetra_vector.at(tetra_id);
        if (it.test(31)) {
            double z = find_polygon_intersection_z(tetra[1], tetra[2], tetra[3]);
            points[i] = z;
            i++;
        }
        if (it.test(30)) {
            double z = find_polygon_intersection_z(tetra[0], tetra[2], tetra[3]);
            points[i] = z;
            i++;
        }
        if (it.test(29)) {
            double z = find_polygon_intersection_z(tetra[0], tetra[1], tetra[3]);
            points[i] = z;
            i++;
        }
        if (it.test(28)) {
            double z = find_polygon_intersection_z(tetra[0], tetra[1], tetra[2]);
            points[i] = z;
            i++;
        }

        if (points.at(0) < points.at(1)) {
            const auto tmp = points.at(0);
            points.at(0) = points.at(1);
            points.at(1) = tmp;
        }

        intersection_data data{.delta_z = points.at(0) - points.at(1),
                               .tetra_id = static_cast<unsigned int>(tetra_id)};

        set_intersection_data set_data{.key_z = points[0], .data = data};

        values.push_back(set_data);
    }

    std::sort(values.begin(), values.end(), set_intersection_data_cmp());

    std::vector<intersection_data> delta_values;
    delta_values.reserve(values.size());

    for (auto& it : values) {
        delta_values.push_back(it.data);
    }

    _intersections_delta = std::move(delta_values);
}

double line::find_polygon_intersection_z(const std::array<double, 3>& p1,
                                         const std::array<double, 3>& p2,
                                         const std::array<double, 3>& p3) {
    /*
     * using matrix plane equation
     *
     * (x - x1) * following minor part
     */
    double x_x1 =
        (_x - p1[0]) * ((p2[1] - p1[1]) * (p3[2] - p1[2]) - (p3[1] - p1[1]) * (p2[2] - p1[2]));

    /*
     * (y - y1) * following minor part
     */
    double y_y1 =
        (_y - p1[1]) * ((p2[0] - p1[0]) * (p3[2] - p1[2]) - (p3[0] - p1[0]) * (p2[2] - p1[2]));

    /*
     * only minor part of (z - z1)
     */
    double z_z1_minor = ((p2[0] - p1[0]) * (p3[1] - p1[1]) - (p3[0] - p1[0]) * (p2[1] - p1[1]));
    double z_value = (y_y1 - x_x1) / z_z1_minor + p1[2];

    return z_value;
}

double line::direct_calculate_ray_value(const std::vector<tetra>& tetra_vector,
                                        tetra_value value_signature) {
    if (_marked_solid) {
        return _mark_value;
    }

    double sum = 0;

    const size_t length = _tetra_intersections.size();
    for (size_t i = 0; i < length; i++) {
        size_t tetra_id = _intersections_delta[i].tetra_id;
        double delta = _intersections_delta[i].delta_z;

        sum = sum + delta * tetra_vector[tetra_id].access_value(value_signature);
    }

    return sum;
}

double line::integrate_ray_value_by_i(const std::vector<tetra>& tetra_vector,
                                      tetra_value alpha_signature, tetra_value q_signature) {
    if (_marked_solid) {
        return _mark_value;
    }

    double I = 0;
    const auto length = static_cast<ssize_t>(_tetra_intersections.size());

    const double limit_alpha = app::instance().config.limit_alpha_value;

    for (ssize_t i = length - 1; i >= 0; i--) {
        size_t tetra_id = _intersections_delta[i].tetra_id;
        double delta = _intersections_delta[i].delta_z;

        double Q = tetra_vector[tetra_id].access_value(q_signature);
        double alpha = tetra_vector[tetra_id].access_value(alpha_signature);

        /*
         * stub to limit huge alpha values raising artifacts
         */
        if (alpha > limit_alpha) {
            alpha = limit_alpha;
        }

        double C = Q - alpha * I;
        if (alpha < std::numeric_limits<double>::epsilon()) {
        } else {
            I = (Q - C * exp(-alpha * delta)) / alpha;
        }
    }
    return I;
}

void line::ts_tetra_intersections_pushback(std::bitset<32> data) {
    std::lock_guard<std::mutex> lock(_tetra_intersections_mutex);
    _tetra_intersections.push_back(data);
}

line::line(line&& ref_val) noexcept {
    _x = ref_val._x;
    _y = ref_val._y;
    _tetra_intersections = std::move(ref_val._tetra_intersections);
    _intersections_delta = std::move(ref_val._intersections_delta);
}

void line::free_memory() {
    std::vector<std::bitset<32>>().swap(_tetra_intersections);
    std::vector<intersection_data>().swap(_intersections_delta);
}

void line::mark_solid_color(double mark_value) {
    _marked_solid = true;
    _mark_value = mark_value;
}