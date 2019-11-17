#include <iostream>
#include <algorithm>
#include <chrono>

#include <lite_tetrahedron.hpp>
#include <line.hpp>
#include <plane.hpp>
#include <application_core.hpp>

constexpr double pi = 3.14159265358979323846;

int main() {
    auto& sys = application_core::instance();

    sys.set_resolution({1900, 1000});
    sys.set_filename("test.vtk");

    std::vector<lite_tetrahedron> tetrahedron_vector{};
    plane current_plane{};
    std::array<double, 4> global_boundaries{};
    try {
        tetrahedron_vector = sys.get_tetrahedron_vector();
        application_core::rotate_tetrahedron_vector(tetrahedron_vector, pi / 12.);
        global_boundaries = application_core::get_global_boundaries(tetrahedron_vector);
        application_core::print_boundaries(global_boundaries);

        current_plane = sys.init_plane_grid(global_boundaries);
    } catch (const std::exception& e) {
        std::cout << "initialization error" << std::endl;
        std::cerr << e.what() << std::endl;
    }

    size_t min_intersections(0), max_intersections(0);
    min_intersections = max_intersections = current_plane.find_intersections_with_tetrahedron(tetrahedron_vector[0], 0);

    auto t1 = std::chrono::high_resolution_clock::now();

    for (size_t i = 1; i < tetrahedron_vector.size(); i++) {
        size_t res = current_plane.find_intersections_with_tetrahedron(tetrahedron_vector[i], i);
        min_intersections = std::min(res, min_intersections);
        max_intersections = std::max(res, max_intersections);
    }

    auto t2 = std::chrono::high_resolution_clock::now();

    auto timer = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    std::cout << "find all intersections complete in: " << timer << " ms" << std::endl;
    std::cout << max_intersections << std::endl;
    std::cout << min_intersections << std::endl;

    std::cout << current_plane.count_all_intersections() << std::endl;

    return 0;
}