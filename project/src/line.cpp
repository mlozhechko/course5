#include <line.hpp>

#include <tuple>

line::line(double x, double y) : _x(x), _y(y) {
    _tetra_intersections.reserve(50);
}

double line::x() {
    return _x;
}

double line::y() {
    return _y;
}

void line::add_tetra_intersection(size_t id, size_t polygon_id) {
    if (!alternation_flag) {
        _tetra_intersections.emplace_back(id);
    }

    _tetra_intersections.back().set(28 + polygon_id);
    alternation_flag = !alternation_flag;
}

size_t line::number_of_intersections() {
    return _tetra_intersections.size();
}

/*
 * 0001 0, 1, 2 points
 * 0010 0, 1, 3 points
 * 0100 0, 2, 3 points
 * 1000 1, 2, 3 points
 */
static std::bitset<32> mask{0x0FFFFFFF};

void line::calculate_intersections(const std::vector<lite_tetrahedron>& tetra_vector) {

    /*
     * {min value, delta value}
     */

    std::vector<std::tuple<double, double, intersection_data>> z_values;
    z_values.reserve(_tetra_intersections.size());

    for (const auto& it: _tetra_intersections) {
//        std::cout << it << std::endl;

        size_t tetra_id = (it & mask).to_ulong();
        std::vector<double> points{};
        points.reserve(2);

        const auto& tetra = tetra_vector.at(tetra_id);
        if (it.test(31)) {
            double z = find_polygon_intersection_z(tetra[1], tetra[2], tetra[3]);
//            std::cout << "1 " << tetra_id << ": (x, y, z) == (" << _x << ", " << _y << ", " << z << ")"<< std::endl;
            points.push_back(z);
        }
        if (it.test(30)) {
            double z = find_polygon_intersection_z(tetra[0], tetra[2], tetra[3]);
//            std::cout << "2 " << tetra_id << ": (x, y, z) == (" << _x << ", " << _y << ", " << z << ")"<< std::endl;
            points.push_back(z);
        }
        if (it.test(29)) {
            double z = find_polygon_intersection_z(tetra[0], tetra[1], tetra[3]);
//            std::cout << "3 " << tetra_id << ": (x, y, z) == (" << _x << ", " << _y << ", " << z << ")"<< std::endl;
            points.push_back(z);
        }
        if (it.test(28)) {
            double z = find_polygon_intersection_z(tetra[0], tetra[1], tetra[2]);
//            std::cout << "4 " << tetra_id << ": (x, y, z) == (" << _x << ", " << _y << ", " << z << ")" << std::endl;
            points.push_back(z);
        }

        if (points.at(0) < points.at(1)) {
            std::swap(points.at(0), points.at(1));
        }

        intersection_data data{
            .delta_z = points.at(0) - points.at(1),
            .tetra_id = static_cast<unsigned int>(tetra_id)
        };
        z_values.emplace_back(points.at(1), points.at(0), data);
    }

    /*
     * TODO:
     * 1. Implement faster custom sort without pair
     */


    std::sort(z_values.begin(), z_values.end(), [](const auto& a, const auto& b) {
        return std::get<0>(a) < std::get<0>(b);
    });

    std::vector<intersection_data> delta_values;
    delta_values.reserve(z_values.size());

    for (auto& it: z_values) {
        delta_values.push_back(std::get<2>(it));

//        std::cout << std::get<0>(it) << " " << std::get<1>(it) << std::endl;
    }
//    std::cout << "--" << std::endl;

    /*
     * TODO:
     * 1. for the god's sake remove debug-tuple from this part of code
     */

    _intersections_delta = std::move(delta_values);
}

double line::find_polygon_intersection_z(const std::array<double, 3>& p1, const std::array<double, 3>& p2,
                                         const std::array<double, 3>& p3) {
    /*
     * using matrix plane equation
     *
     * (x - x1) * following minor part
     */
    double x_x1 = (_x - p1[0]) * ( (p2[1] - p1[1]) * (p3[2] - p1[2]) - (p3[1] - p1[1]) * (p2[2] - p1[2]) );

    /*
     * (y - y1) * following minor part
     */
    double y_y1 = (_y - p1[1]) * ( (p2[0] - p1[0]) * (p3[2] - p1[2]) - (p3[0] - p1[0]) * (p2[2] - p1[2]) );

    /*
     * only minor part of (z - z1)
     */
    double z_z1_minor = ( (p2[0] - p1[0]) * (p3[1] - p1[1]) - (p3[0] - p1[0]) * (p2[1] - p1[1]) );
    double z_value = (y_y1 - x_x1) / z_z1_minor + p1[2];

    return z_value;
}

double line::direct_calculate_ray_value(const std::vector<lite_tetrahedron>& tetra_vector, tetra_value value_signature) {
    double sum = 0;

    size_t length = _tetra_intersections.size();
    for (size_t i = 0; i < length; i++) {
        size_t tetra_id = _intersections_delta.at(i).tetra_id;
        double delta = _intersections_delta.at(i).delta_z;

        sum = sum + delta * tetra_vector.at(tetra_id).access_value(value_signature);
    }

    return sum;
}

double line::integrate_ray_value_by_i(const std::vector<lite_tetrahedron>& tetra_vector, tetra_value alpha_signature,
                                      tetra_value q_signature) {
    double I = 0;

    size_t length = static_cast<ssize_t>(_tetra_intersections.size());
    for (ssize_t i = length - 1; i >= 0; i--) {
        size_t tetra_id = _intersections_delta.at(i).tetra_id;
        double delta = _intersections_delta.at(i).delta_z;

        double Q = tetra_vector.at(tetra_id).access_value(q_signature);
        double alpha = tetra_vector.at(tetra_id).access_value(alpha_signature);

        double C = Q - alpha * I;
        if (alpha < std::numeric_limits<double>::epsilon()) {
//            std::cout << "warning";
        } else {
            I = (Q - C * exp(-alpha * delta)) / alpha;
        }
//        std::cout << I << " ";
    }
//    std::cout << std::endl;
    return I;
}
