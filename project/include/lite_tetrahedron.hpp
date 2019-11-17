#pragma once

#include <array>
#include <cmath>
#include <string>
#include <iostream>

class lite_tetrahedron {
public:
    explicit lite_tetrahedron(const std::array<std::array<double, 3>, 4>& points, double alpha, double q);

    const std::array<double, 3>& operator[](size_t i) const {
        return _points[i];
    }

    /*
     * rotation by angle by y axis
     */
    void rotate_x(double angle);

    /*
     * get tetrahedrons boundaries
     * std::array<double, 4> = {x_max, x_min, y_max, y_min}
     */
    std::array<double, 4> get_boundaries();

private:
    static void point_rotate_x(double alpha, std::array<double, 3>& point);
    friend std::ostream& operator<<(std::ostream& os, const lite_tetrahedron& lt);

    std::array<std::array<double, 3>, 4> _points{};
    double _q;
    double _alpha;
};

std::ostream& operator<<(std::ostream& os, const lite_tetrahedron& lt);