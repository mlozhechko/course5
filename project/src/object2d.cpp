#include <object2d.hpp>

object2d::object2d(const std::pair<double_matrix, double_matrix>& object2d_data) {
    _object2d_data = object2d_data;
}

void object2d::export_to_vti(const std::string& filename) {
    const size_t res_x = _object2d_data.first.size();
    const size_t res_y = _object2d_data.first.at(0).size();

    vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
    imageData->SetDimensions(res_x, res_y, 1);
    imageData->AllocateScalars(VTK_DOUBLE, 2);

    int* dims = imageData->GetDimensions();

    for (int y = 0; y < dims[1]; y++) {
        for (int x = 0; x < dims[0]; x++) {
            auto* pixel = static_cast<double*>(imageData->GetScalarPointer(x, y, 0));
            pixel[0] = _object2d_data.first[x][y];
            pixel[1] = _object2d_data.second[x][y];
        }
    }

    vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(imageData);
    writer->Write();
}
