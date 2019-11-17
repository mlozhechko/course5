#include <line.hpp>


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

    _tetra_intersections.back().set(31 - polygon_id);
    alternation_flag = !alternation_flag;
}

size_t line::number_of_intersections() {
    return _tetra_intersections.size();
}