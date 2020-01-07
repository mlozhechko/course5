#pragma once

#include <array>
#include <cmath>
#include <string>
#include <iostream>

enum class tetra_value : size_t {
    alpha = 0, Q = 1
};

class tetra {
public:
    explicit tetra(const std::array<std::array<double, 3>, 4>& points, double alpha, double q);

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

    double access_value(tetra_value value) const;

    /* debug */
    double delta_z();

private:
    static void point_rotate_x(double alpha, std::array<double, 3>& point);
    friend std::ostream& operator<<(std::ostream& os, const tetra& lt);

    std::array<std::array<double, 3>, 4> _points{};
    std::array<double, 2> _tetra_values;
};

std::ostream& operator<<(std::ostream& os, const tetra& lt);