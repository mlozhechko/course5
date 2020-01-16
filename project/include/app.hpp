#pragma once

#include <array>
#include <atomic>
#include <chrono>
#include <string>
#include <thread>
#include <vector>

#include <config.hpp>
#include <plane.hpp>
#include <tetra.hpp>

#include <omp.h>

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkImageData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkXMLImageDataWriter.h>

namespace app {
vtkSmartPointer<vtkUnstructuredGrid> init_vtk_grid(const std::string& filename);
std::vector<tetra> get_tetrahedron_vector(const std::string& filename);
std::vector<tetra> get_tetrahedron_vector_binary(const std::string& filename);
void rotate_tetrahedron_vector(std::vector<tetra>&, double angle);
std::array<double, 4> get_domain_boundaries(std::vector<tetra>&);
void print_boundaries(const std::array<double, 4>&);
//plane init_plane_grid(size_t res_x, size_t res_y, std::array<double, 4>& global_boundaries);

void find_tetrahedron_vector_intersections_with_lines(const std::vector<tetra>& tetrahedron_vector,
                                                      plane& task_plane);

std::pair<float_matrix, float_matrix> trace_rays(plane& current_plane,
                                                 std::vector<tetra>& tetrahedron_vector,
                                                 tetra_value value_alpha, tetra_value value_q);

void produce_result(std::pair<float_matrix, float_matrix>& data, size_t res_x, size_t res_y,
                    const std::string& result_filename);
} // namespace app
