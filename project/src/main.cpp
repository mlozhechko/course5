#include <iostream>
#include <algorithm>

#include <config.hpp>
#include <tetra.hpp>
#include <line.hpp>
#include <plane.hpp>
#include <app.hpp>


int main() {
    std::cout << "[1] ";
    std::string filename = "test2.vtk";
    std::string res_filename = "result5.vti";
    size_t res_x = 1000;
    size_t res_y = 1000;
    double projection_angle = 5;


    std::cout << "defined grid resolution: " << res_x << "x" << res_y << std::endl;
    std::cout << "source file: " << filename << std::endl;

    std::vector<tetra> tetrahedron_vector{};
    plane current_plane{};
    std::array<double, 4> domain_boundaries{};
    try {
        auto t1 = std::chrono::high_resolution_clock::now();

        tetrahedron_vector = app::get_tetrahedron_vector_binary(filename);
        app::rotate_tetrahedron_vector(tetrahedron_vector, pi * projection_angle / 180.);
        domain_boundaries = app::get_domain_boundaries(tetrahedron_vector);
        current_plane = app::init_plane_grid(res_x, res_y, domain_boundaries);

        auto t2 = std::chrono::high_resolution_clock::now();

        auto timer = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        std::cout << "initialization complete in: " << timer << " ms" << std::endl;
    } catch (const std::exception& e) {
        std::cerr
            << "initialization error" << std::endl
            << e.what() << std::endl;

        return -1;
    }
    std::cout << "data loaded from file, tetrahedrons rotated. ready for ray-tracing" << std::endl;

    std::cout << "[2] ";
    try {
        app::find_tetrahedron_vector_intersections_with_lines(tetrahedron_vector, current_plane);
    } catch (const std::exception& e) {
        std::cerr
            << "iteration through tetrahedron array failure" << std::endl
            << e.what() << std::endl;
        return -2;
    }
    std::cout << "rays and tetrahedrons matching completed" << std::endl;
    std::cout << "[3] ";

    std::pair<float_matrix, float_matrix> result;
    try {
        result = app::trace_rays(current_plane, tetrahedron_vector, tetra_value::alpha, tetra_value::Q);
    } catch (const std::exception& e) {
        std::cerr
            << "ray tracing failure" << std::endl
            << e.what() << std::endl;
        return -3;
    }
    std::cout << "ray tracing complete" << std::endl;

    std::cout << "[4] ";
    try {
        app::produce_result(result, res_x, res_y, res_filename);
    } catch (const std::exception& e) {
        std::cerr
            << "producing result failure" << std::endl
            << e.what() << std::endl;
        return -4;
    }
    std::cout << "result file: " << res_filename << std::endl;
    std::cout << "producing result representation completed" << std::endl;

    return 0;
}