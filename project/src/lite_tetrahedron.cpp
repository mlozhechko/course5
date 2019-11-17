#include <lite_tetrahedron.hpp>

lite_tetrahedron::lite_tetrahedron(const std::array<std::array<double, 3>, 4>& points, double alpha, double q)
: _points(points), _q(q), _alpha(alpha) {}

void lite_tetrahedron::rotate_x(double angle) {
    for (auto& i: _points) {
        point_rotate_x(angle, i);
    }
}

std::array<double, 4> lite_tetrahedron::get_boundaries() {
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

        flag = static_cast<double>(res[1] >_points[i][0]);
        res[1] = flag * _points[i][0] + (1 - flag) * res[1];

        flag = static_cast<double>(res[2] < _points[i][1]);
        res[2] = flag * _points[i][1] + (1 - flag) * res[2];

        flag = static_cast<double>(res[3] > _points[i][1]);
        res[3] = flag * _points[i][1] + (1 - flag) * res[3];
    }

    return res;
}

void lite_tetrahedron::point_rotate_x(double alpha, std::array<double, 3>& point) {
    double tmp_point1 = point[1];
    point[1] = point[1] * cos(alpha) - point[2] * sin(alpha);
    point[2] = tmp_point1 * sin(alpha) + point[2] * cos(alpha);
}

std::ostream& operator<<(std::ostream& os, const lite_tetrahedron& lt)
{
    os << "lite tetrahedron. q = " << lt._q << ", alpha = " << lt._alpha << std::endl;
    for (const auto& it: lt._points) {
        os << it[0] << " " << it[1] << " " << it[2] << std::endl;
    }
    return os;
}