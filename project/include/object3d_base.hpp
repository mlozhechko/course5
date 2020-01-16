#pragma once

#include <memory>
#include <vector>
#include <functional>
#include <array>

#include <tetra.hpp>
#include <config.hpp>

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkImageData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>

/*
 * 3d object consists of tetra
 */
class object3d_base {
public:
    object3d_base() = default;
    virtual ~object3d_base() = default;

    void read_vtk_file(const std::string& filename, const std::vector<std::string>& scalar_labels);
    //virtual void write_to_vtk_file();
    void init_polar(const std::function<double(std::array<double, 3>)>& potential_function, double x0, double y0,
                    double z0, double level_value, double step, double angle_step);

    std::shared_ptr<std::vector<tetra>> get_pointer();

    virtual void rotate_around_x_axis(double angle);
    virtual std::array<double, 4> get_boundaries();

protected:
    std::shared_ptr<std::vector<tetra>> _data = std::make_shared<std::vector<tetra>>();
};