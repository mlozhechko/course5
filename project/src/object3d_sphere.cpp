#include <object3d_sphere.hpp>

static double vector_2_norm(const std::array<double, 3>& vec) {
    double sum{};
    for (const auto& i : vec) {
        sum += i * i;
    }
    return sqrt(sum);
}

object3d_sphere::object3d_sphere(const point& center, double R) {
    auto sphere_func = [&](const point& p){
        point internal_point{p[0] - center[0], p[1] - center[1], p[2] - center[2]};
        return vector_2_norm(internal_point);
    };

    init_polar(sphere_func, center[0], center[1], center[2], R, 0.001, 256, tetra_type::solid, ACC_DISK_COLOR);
}