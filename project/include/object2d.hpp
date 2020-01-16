#pragma once

#include <utility>
#include <vector>
#include <string>

#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkXMLImageDataWriter.h>

using float_matrix = std::vector<std::vector<float>>;

class object2d {
public:
    explicit object2d(const std::pair<float_matrix, float_matrix>& object2d_data);
    void export_to_vti(const std::string& filename);

    /*
     * class can be filled with exports to to other data formats
     */
private:
    std::pair<float_matrix, float_matrix> _object2d_data;
};