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

struct config_str {
    std::string file = {};
    std::string destination = {};
    int threads = 1;
    size_t resolution_x = 1100;
    size_t resolution_y = 900;
} config;

namespace po = boost::program_options;
bool program_options(int argc, char **argv) {
    bool is_ready = false;

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message")
        ("file,f", po::value<std::string>(), "source file")
        ("destination,d", po::value<std::string>(),"destination file")
        ("threads,j", po::value<int>(), "amount of parallel threads")
        ("resolution_x,x", po::value<size_t>(), "set x axis resolution")
        ("resolution_y,y", po::value<size_t>(), "set y axis resolution");

    po::variables_map vm;
    store(parse_command_line(argc, argv, desc), vm);
    notify(vm);

    if (vm.count("help")) {
        std::cout << desc << std::endl;
        return is_ready;
    }

    if (vm.count("file") && vm.count("destination")) {
        config.file = vm["file"].as<std::string>();
        config.destination = vm["destination"].as<std::string>();
    } else {
        std::cout << "Error! Source filename and destination filename must be specified" << std::endl;
        std::cout << desc << std::endl;
        return is_ready;
    }
    is_ready = true;

    if (vm.count("threads")) {
        config.threads = vm["threads"].as<int>();
    } else {
        config.threads = std::max(static_cast<int>(std::thread::hardware_concurrency()), 1);
    }

    if (vm.count("resolution_x")) {
        config.resolution_x = vm["resolution_x"].as<size_t>();
    }

    if (vm.count("resolution_y")) {
        config.resolution_y = vm["resolution_y"].as<size_t>();
    }

    return is_ready;
}

int main(int argc, char** argv) {
   bool is_ready = program_options(argc, argv);
   if (!is_ready) {
       return 0;
   }

    /*
     * {x_max, x_min, y_max, y_min}
     */
    std::vector<double> domain = {2.2, -0.2, 0.9, -0.9};
    double projection_angle = 0;

    std::cout << "defined grid resolution: " << config.resolution_x << "x" << config.resolution_y << std::endl;
    std::cout << "source file: " << config.file << std::endl;
    std::cout << "amount of parallel threads: " << config.threads << std::endl;

    auto acc_disk = object3d_accretion_disk{config.file};
    acc_disk.rotate_around_x_axis(PI * projection_angle / 180.);

    auto roche_lobe = object3d_roche_lobe{{ACC_X0, ACC_Y0, ACC_Z0}, L, 0, M_ACC, M_DONOR, OMEGA};
    auto acc_sphere = object3d_sphere{{ACC_X0, ACC_Y0, ACC_Z0}, ACC_DISK_R};

    plane base_plane{config.resolution_x, config.resolution_y, {acc_disk, roche_lobe, acc_sphere}, domain};
    base_plane.find_intersections();
    object2d result = base_plane.trace_rays(tetra_value::alpha, tetra_value::Q);
    result.export_to_vti(config.destination);

    return 0;
}