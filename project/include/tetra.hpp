#pragma once

#include <array>
#include <cmath>
#include <iostream>
#include <string>

enum class tetra_value : size_t { alpha = 0, solid_color = 0, Q = 1 };

enum class tetra_type {transparent = 0, solid = 1};

class tetra {
public:
    explicit tetra(const std::array<std::array<double, 3>, 4>& points, double v1, double v2, tetra_type ttype = tetra_type::transparent);

    const std::array<double, 3>& operator[](size_t i) const {
        return _points[i];
    }

    void rotate_around_x_axis(double angle);
    void rotate_around_y_axis(double angle, double x0);

    /*
     * get tetrahedrons boundaries
     * std::array<double, 4> = {x_max, x_min, y_max, y_min}
     */
    std::array<double, 4> get_boundaries();

    double access_value(tetra_value value) const;

    tetra_type get_tetra_type() const;

    /* debug */
    double delta_z();

private:
    static void point_rotate_around_x_axis(double alpha, std::array<double, 3>& point);
    static void point_rotate_around_y_axis(double alpha, std::array<double, 3>& point, double x0 = 0);

    friend std::ostream& operator<<(std::ostream& os, const tetra& lt);

    std::array<std::array<double, 3>, 4> _points{};
    std::array<double, 2> _tetra_values;

    tetra_type _tetra_type;
};

std::ostream& operator<<(std::ostream& os, const tetra& lt);