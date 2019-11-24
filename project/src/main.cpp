#include <iostream>
#include <algorithm>
#include <chrono>

#include <lite_tetrahedron.hpp>
#include <line.hpp>
#include <plane.hpp>
#include <application_core.hpp>
#include <vtkImageData.h>
#include <vtkXMLImageDataWriter.h>

constexpr double pi = 3.14159265358979323846;

int main() {
    auto& sys = application_core::instance();

    size_t x_res = 2000;
    size_t y_res = 1000;

    sys.set_resolution({x_res, y_res});
    sys.set_filename("test.vtk");

    std::vector<lite_tetrahedron> tetrahedron_vector{};
    plane current_plane{};
    std::array<double, 4> global_boundaries{};
    try {
        tetrahedron_vector = sys.get_tetrahedron_vector();
        application_core::rotate_tetrahedron_vector(tetrahedron_vector, pi * 85. / 180.);
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

//    size_t res = current_plane.find_intersections_with_tetrahedron(tetrahedron_vector[610], 610);

//    current_plane.print_all_lines_with_intersection();

    auto t2 = std::chrono::high_resolution_clock::now();

    auto timer = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    std::cout << "find all intersections complete in: " << timer << " ms" << std::endl;
    std::cout << max_intersections << std::endl;
    std::cout << min_intersections << std::endl;

    t1 = std::chrono::high_resolution_clock::now();
    std::cout << current_plane.count_all_intersections() << std::endl;

    std::vector<std::vector<float>> result_alpha = current_plane.trace_rays(tetrahedron_vector);

    t2 = std::chrono::high_resolution_clock::now();
    timer = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();

    std::cout << "trace all rays complete in: " << timer << " ms" << std::endl;

    vtkSmartPointer<vtkImageData> imageData =
        vtkSmartPointer<vtkImageData>::New();
    imageData->SetDimensions(x_res, y_res, 1);
    imageData->AllocateScalars(VTK_DOUBLE, 2);

    int* dims = imageData->GetDimensions();

    // Fill every entry of the image data with "2.0"
    for (int z = 0; z < dims[2]; z++)
    {
        for (int y = 0; y < dims[1]; y++)
        {
            for (int x = 0; x < dims[0]; x++)
            {
                double* pixel = static_cast<double*>(imageData->GetScalarPointer(x,y,z));
                pixel[0] = result_alpha[x][y];
                pixel[1] = 1.0;
            }
        }
    }

    vtkSmartPointer<vtkXMLImageDataWriter> writer =
        vtkSmartPointer<vtkXMLImageDataWriter>::New();
    writer->SetFileName("testfile.vti");
    writer->SetInputData(imageData);
    writer->Write();

    return 0;
}