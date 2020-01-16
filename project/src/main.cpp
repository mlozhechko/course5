#include <algorithm>
#include <iostream>

#include <object3d_accretion_disk.hpp>
#include <object3d_roche_lobe.hpp>

#include <app.hpp>
#include <config.hpp>
#include <line.hpp>
#include <plane.hpp>
#include <tetra.hpp>
#include <object2d.hpp>

int main() {
    std::string filename = "test2.vtk";
    std::string res_filename = "result5.vti";
    size_t res_x = 1400;
    size_t res_y = 1000;
    double projection_angle = -5;

    std::cout << "defined grid resolution: " << res_x << "x" << res_y << std::endl;
    std::cout << "source file: " << filename << std::endl;

    auto acc_disk = object3d_accretion_disk{filename};
    acc_disk.rotate_around_x_axis(PI * projection_angle / 180.);

    auto roche_lobe = object3d_roche_lobe{{acc_x0, acc_y0, acc_z0}, L, 0, Ma, Md, OMEGA};

    plane base_plane{res_x, res_y, {acc_disk, roche_lobe}};
    base_plane.find_intersections();
    object2d result = base_plane.trace_rays(tetra_value::alpha, tetra_value::Q);
    result.export_to_vti(res_filename);

//    std::vector<tetra> tetrahedron_vector{};
//    std::vector<tetra> roche_lobe{};
//    plane current_plane{};
//    std::array<double, 4> domain_boundaries{};
//    try {
//        auto t1 = std::chrono::high_resolution_clock::now();
//        tetrahedron_vector = app::get_tetrahedron_vector_binary(filename);
//        auto t2 = std::chrono::high_resolution_clock::now();
//
//        app::rotate_tetrahedron_vector(tetrahedron_vector, PI * projection_angle / 180.);
//
//        domain_boundaries = app::get_domain_boundaries(tetrahedron_vector);
//        //        std::cout << domain_boundaries[0] << " " << domain_boundaries[1] << " " <<
//        //        domain_boundaries[2] << " " << domain_boundaries[3] << std::endl;
//        //        domain_boundaries = {1.692, 0.308021, 0.689701, -0.689681};
//        domain_boundaries[1] -= 0.5;
//        domain_boundaries[0] += 0.5;
//
//        current_plane = app::init_plane_grid(res_x, res_y, domain_boundaries);
//        roche_lobe = current_plane.build3d_model_donor_roche_lobe();
//
//        std::cout << "roche lobe tetra size: " << roche_lobe.size() << std::endl;
//        std::cout << "tetra vector size: " << tetrahedron_vector.size() << std::endl;
//        app::rotate_tetrahedron_vector(roche_lobe, PI * projection_angle / 180.);
//
//        auto timer = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
//        std::cout << "initialization complete in: " << timer << " ms" << std::endl;
//    } catch (const std::exception& e) {
//        std::cerr << "initialization error" << std::endl << e.what() << std::endl;
//
//        return -1;
//    }
//    std::cout << "data loaded from file, tetrahedrons rotated. ready for ray-tracing" << std::endl;
//
//    std::cout << "[2] ";
//    try {
//        tetrahedron_vector.insert(tetrahedron_vector.end(), roche_lobe.begin(), roche_lobe.end());
//        app::find_tetrahedron_vector_intersections_with_lines(tetrahedron_vector, current_plane);
//    } catch (const std::exception& e) {
//        std::cerr << "iteration through tetrahedron array failure" << std::endl
//                  << e.what() << std::endl;
//        return -2;
//    }
//    std::cout << "rays and tetrahedrons matching completed" << std::endl;
//    std::cout << "[3] ";
//
//    std::pair<float_matrix, float_matrix> result;
//    try {
//        result =
//            app::trace_rays(current_plane, tetrahedron_vector, tetra_value::alpha, tetra_value::Q);
//    } catch (const std::exception& e) {
//        std::cerr << "ray tracing failure" << std::endl << e.what() << std::endl;
//        return -3;
//    }
//    std::cout << "ray tracing complete" << std::endl;
//
//    std::cout << "[4] ";
//    try {
//        app::produce_result(result, res_x, res_y, res_filename);
//    } catch (const std::exception& e) {
//        std::cerr << "producing result failure" << std::endl << e.what() << std::endl;
//        return -4;
//    }
//    std::cout << "result file: " << res_filename << std::endl;
//    std::cout << "producing result representation completed" << std::endl;

    return 0;
}