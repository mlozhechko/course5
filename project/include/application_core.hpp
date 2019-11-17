#pragma once

#include <string>
#include <vector>
#include <array>

#include <plane.hpp>
#include <lite_tetrahedron.hpp>

#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>

/*
 * application singleton class
 */
class application_core{
public:
    static application_core& instance();
    std::vector<lite_tetrahedron> get_tetrahedron_vector();
    static void rotate_tetrahedron_vector(std::vector<lite_tetrahedron>&, double angle);
    static std::array<double, 4> get_global_boundaries(std::vector<lite_tetrahedron>&);
    static void print_boundaries(const std::array<double, 4>&);
    plane init_plane_grid(std::array<double, 4>&);

    application_core(application_core&) = delete;
    application_core(application_core&&) = delete;
    void operator=(application_core&) = delete;
    void operator=(application_core&&) = delete;

    void set_filename(const std::string&);
    void set_resolution(const std::array<size_t, 2>&);

private:
    std::string _filename{};
    size_t _res_x{0}, _res_y{0};

    application_core();
    static vtkSmartPointer<vtkUnstructuredGrid> init_vtk_grid(const std::string& filename);

};


