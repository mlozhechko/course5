#include <iostream>
#include <algorithm>

#include <lite_tetrahedron.hpp>
#include <line.hpp>
#include <plane.hpp>
#include <app.hpp>

constexpr double pi = 3.14159265358979323846;

int main() {
    std::string filename = "test.vtk";
    std::string res_filename = "result61.vti";
    size_t res_x = 200;
    size_t res_y = 200;
    double projection_angle = 61;

    std::cout << "defined grid resolution: " << res_x << "x" << res_y << std::endl;
    std::cout << "source file: " << filename << std::endl;

    std::vector<lite_tetrahedron> tetrahedron_vector{};
    plane current_plane{};
    std::array<double, 4> domain_boundaries{};
    try {
        tetrahedron_vector = app::get_tetrahedron_vector(filename);
        app::rotate_tetrahedron_vector(tetrahedron_vector, pi * projection_angle / 180.);
        domain_boundaries = app::get_domain_boundaries(tetrahedron_vector);
        current_plane = app::init_plane_grid(res_x, res_y, domain_boundaries);
    } catch (const std::exception& e) {
        std::cerr
            << "initialization error" << std::endl
            << e.what() << std::endl;

        return -1;
    }
    std::cout << "[1] grid initialized, source data produced" << std::endl;

    try {
        app::find_tetrahedron_vector_intersections_with_lines(tetrahedron_vector, current_plane);
    } catch (const std::exception& e) {
        std::cerr
            << "iteration through tetrahedron array failure" << std::endl
            << e.what() << std::endl;
        return -2;
    }
    std::cout << "[2] rays and tetrahedrons matching completed" << std::endl;

    std::vector<std::vector<float>> result_alpha;
    try {
        result_alpha = app::direct_trace_rays(current_plane, tetrahedron_vector, tetra_value::alpha);
    } catch (const std::exception& e) {
        std::cerr
            << "ray tracing failure" << std::endl
            << e.what() << std::endl;
        return -3;
    }
    std::cout << "[3] ray tracing complete" << std::endl;

    try {
        app::produce_result(result_alpha, res_filename, res_x, res_y);
    } catch (const std::exception& e) {
        std::cerr
            << "producing result failure" << std::endl
            << e.what() << std::endl;
        return -4;
    }
    std::cout << "result file: " << res_filename << std::endl;
    std::cout << "[4] producing result representation completed" << std::endl;

    return 0;
}