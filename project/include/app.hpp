#pragma once

#include <string>
#include <vector>
#include <array>
#include <chrono>

#include <plane.hpp>
#include <lite_tetrahedron.hpp>

#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkImageData.h>
#include <vtkXMLImageDataWriter.h>

namespace app {
    vtkSmartPointer<vtkUnstructuredGrid> init_vtk_grid(const std::string& filename);
    std::vector<lite_tetrahedron> get_tetrahedron_vector(const std::string& filename);
    void rotate_tetrahedron_vector(std::vector<lite_tetrahedron>&, double angle);
    std::array<double, 4> get_domain_boundaries(std::vector<lite_tetrahedron>&);
    void print_boundaries(const std::array<double, 4>&);
    plane init_plane_grid(size_t res_x, size_t res_y, std::array<double, 4>& global_boundaries);

    void find_tetrahedron_vector_intersections_with_lines(
        const std::vector<lite_tetrahedron>& tetrahedron_vector,
        plane& task_plane);

    std::vector<std::vector<float>> direct_trace_rays(
        plane& current_plane,
        std::vector<lite_tetrahedron>& tetrahedron_vector,
        tetra_value value);

    void
    produce_result(std::vector<std::vector<float>>& data_x, const std::string& result_filename, size_t res_x, size_t res_y);
}


