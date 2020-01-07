#include <tetra.hpp>

tetra::tetra(const std::array<std::array<double, 3>, 4>& points, double alpha, double q)
    : _points(points) {
    _tetra_values[static_cast<size_t>(tetra_value::alpha)] = alpha;
    _tetra_values[static_cast<size_t>(tetra_value::Q)] = q;
}

void tetra::rotate_x(double angle) {
    for (auto& i: _points) {
        point_rotate_x(angle, i);
    }
}

std::array<double, 4> tetra::get_boundaries() {
    /*
     * TODO:
     * 1. remove(!)
     */
    std::array<double, 4> res{};
    res[1] = res[0] = _points[0][0];
    res[3] = res[2] = _points[0][1];

    for (size_t i = 1; i < 4; i++) {
        /*
         * GPU optimization to avoid if statements (blank)
         *
         */
        auto flag = static_cast<double>(res[0] < _points[i][0]);
        res[0] = flag * _points[i][0] + (1 - flag) * res[0];

        flag = static_cast<double>(res[1] > _points[i][0]);
        res[1] = flag * _points[i][0] + (1 - flag) * res[1];

        flag = static_cast<double>(res[2] < _points[i][1]);
        res[2] = flag * _points[i][1] + (1 - flag) * res[2];

        flag = static_cast<double>(res[3] > _points[i][1]);
        res[3] = flag * _points[i][1] + (1 - flag) * res[3];
    }

    return res;
}

void tetra::point_rotate_x(double alpha, std::array<double, 3>& point) {
    double tmp_point1 = point[1];
    point[1] = point[1] * cos(alpha) - point[2] * sin(alpha);
    point[2] = tmp_point1 * sin(alpha) + point[2] * cos(alpha);
}

std::ostream& operator<<(std::ostream& os, const tetra& lt) {
    double alpha = lt.access_value(tetra_value::alpha);
    double q = lt.access_value(tetra_value::Q);

    os << "lite tetrahedron. q = " << q << ", alpha = " << alpha << std::endl;
    for (const auto& it: lt._points) {
        os << it[0] << " " << it[1] << " " << it[2] << std::endl;
    }
    return os;
}

double tetra::delta_z() {
    double max_z = _points[0][2];
    double min_z = _points[0][2];
    for (size_t i = 1; i < 4; i++) {
        max_z = std::max(_points[i][2], max_z);
        min_z = std::min(_points[i][2], min_z);
    }
    return max_z - min_z;
}

double tetra::access_value(tetra_value value) const {
    return _tetra_values[static_cast<size_t>(value)];
}
