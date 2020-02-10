#pragma once

#include <string>
#include <utility>
#include <vector>

#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLImageDataWriter.h>

using double_matrix = std::vector<std::vector<double>>;

class object2d {
public:
    explicit object2d(const std::pair<double_matrix, double_matrix>& object2d_data);
    void export_to_vti(const std::string& filename);

    /*
     * class can be filled with exports to to other data formats
     */
private:
    std::pair<double_matrix, double_matrix> _object2d_data;
};