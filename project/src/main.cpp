#include <algorithm>
#include <iostream>

#include <object3d_accretion_disk.hpp>
#include <object3d_roche_lobe.hpp>
#include <object3d_sphere.hpp>

#include <config.hpp>
#include <line.hpp>
#include <plane.hpp>
#include <tetra.hpp>
#include <object2d.hpp>

#include <boost/program_options.hpp>

namespace po = boost::program_options;
void program_options() {
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message")
        ("resolution_x,rx", po::value<int>(), "set x axis resolution")
        ("resolution_y,ry", po::value<int>(), "set y axis resolution");
}

int main(int argc, char** argv) {
    program_options();

    std::string filename = "test2.vtk";
    std::string res_filename = "result5.vti";

    size_t res_x = 2200;
    size_t res_y = 1800;
    /*
     * {x_max, x_min, y_max, y_min}
     */
    std::vector<double> domain = {2.2, -0.2, 0.9, -0.9};
    double projection_angle = 0;

    std::cout << "defined grid resolution: " << res_x << "x" << res_y << std::endl;
    std::cout << "source file: " << filename << std::endl;

    auto acc_disk = object3d_accretion_disk{filename};
    acc_disk.rotate_around_x_axis(PI * projection_angle / 180.);

    auto roche_lobe = object3d_roche_lobe{{ACC_X0, ACC_Y0, ACC_Z0}, L, 0, M_ACC, M_DONOR, OMEGA};
    auto acc_sphere = object3d_sphere{{ACC_X0, ACC_Y0, ACC_Z0}, ACC_DISK_R};

    plane base_plane{res_x, res_y, {acc_disk, roche_lobe, acc_sphere}, domain};
    base_plane.find_intersections();
    object2d result = base_plane.trace_rays(tetra_value::alpha, tetra_value::Q);
    result.export_to_vti(res_filename);

    return 0;
}