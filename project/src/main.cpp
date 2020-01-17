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
bool program_options(int argc, char **argv) {
    auto& config = app::instance().config;
    bool is_ready = false;

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message")
        ("file,f", po::value<std::string>(), "source file")
        ("destination,d", po::value<std::string>(),"destination file")
        ("threads,j", po::value<int>(), "number of parallel threads")
        ("resolution_x,x", po::value<size_t>()->default_value(1200), "set x axis resolution")
        ("resolution_y,y", po::value<size_t>()->default_value(900), "set y axis resolution")
        ("angle_around_x,X", po::value<double>()->default_value(0.), "rotate view plane by angle around x axis")
        ("angle_around_y,Y", po::value<double>()->default_value(0.), "rotate view plane by angle around y axis")
        ("donor_angle,D", po::value<double>()->default_value(0.), "initial donor angle around y axis");

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

    config.resolution_x = vm["resolution_x"].as<size_t>();
    config.resolution_y = vm["resolution_y"].as<size_t>();
    config.angle_around_x = vm["angle_around_x"].as<double>();
    config.angle_around_y = vm["angle_around_y"].as<double>();
    config.donor_angle = vm["donor_angle"].as<double>();

    return is_ready;
}

const auto& timestamp = std::chrono::high_resolution_clock::now;

int main(int argc, char** argv) {
    const auto& config = app::instance().config;
    bool is_ready = program_options(argc, argv);
    if (!is_ready) {
       return 0;
    }

    /*
     * {x_max, x_min, y_max, y_min}
     */
    std::vector<double> domain = {2.2, -0.2, 0.9, -0.9};

    std::cout << "Defined grid resolution: " << config.resolution_x << "x" << config.resolution_y << std::endl;
    std::cout << "Source file: " << config.file << std::endl;
    std::cout << "Number of parallel threads: " << config.threads << std::endl;
    std::cout << "Initial rotate angle of roche lobe: " << config.donor_angle << " Pi" << std::endl;
    std::cout << "Plane angle around x: " << config.angle_around_x << " Pi" << std::endl;
    std::cout << "Plane angle around y: " << config.angle_around_y << " Pi" << std::endl;
    omp_set_num_threads(config.threads / 2);

    auto t1 = timestamp();
    object3d_accretion_disk acc_disk{};
    std::thread acc_t([&]() {
        /*
         * VTK library is very slow in terms of working with large file
         * therefore it's better to put this work in parallel thread
         */

        acc_disk = object3d_accretion_disk{config.file};
        acc_disk.rotate_around_x_axis(config.angle_around_x * PI);
        acc_disk.rotate_around_y_axis(config.angle_around_y * PI, ACC_X0);
    });

    auto roche_lobe = object3d_roche_lobe{{ACC_X0, ACC_Y0, ACC_Z0}, L, config.donor_angle * PI, M_ACC, M_DONOR, OMEGA};
    roche_lobe.rotate_around_x_axis(config.angle_around_x * PI);
    roche_lobe.rotate_around_y_axis(config.angle_around_y * PI, ACC_X0);

    auto acc_sphere = object3d_sphere{{ACC_X0, ACC_Y0, ACC_Z0}, ACC_DISK_R};
    acc_t.join();
    auto t2 = timestamp();
    std::cout
        << "Loading data with VTK lib and other preparations completed in "
        << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count()
        << " ms. "
        << std::endl;

    omp_set_num_threads(config.threads);
    t1 = timestamp();
    plane base_plane{config.resolution_x, config.resolution_y, {acc_disk, roche_lobe, acc_sphere}, domain};
    base_plane.find_intersections();
    object2d result = base_plane.trace_rays(tetra_value::alpha, tetra_value::Q);
    t2 = timestamp();
    std::cout
        << "Ray-tracing completed in "
        << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count()
        << " ms. "
        << std::endl;

    result.export_to_vti(config.destination);
    std::cout << "Result exported. Calculations completed." << std::endl;

    return 0;
}