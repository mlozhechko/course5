#include <app.hpp>


vtkSmartPointer<vtkUnstructuredGrid> app::init_vtk_grid(const std::string& filename) {
    auto reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader->SetFileName(filename.c_str());
    reader->SetReadAllScalars(true);
    reader->Update();

//    reader->Print(std::cout);
    return reader->GetOutput();
}

std::vector<tetra> app::get_tetrahedron_vector(const std::string& filename) {
    std::vector<tetra> tetrahedron_vector{};
    vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid{init_vtk_grid(filename)};
    vtkSmartPointer<vtkCellData> scalar_data = unstructured_grid->GetCellData();

    size_t number_of_cells = unstructured_grid->GetNumberOfCells();

    vtkSmartPointer<vtkDataArray> scalars_alpha = scalar_data->GetScalars("alpha");
    vtkSmartPointer<vtkDataArray> scalars_q = scalar_data->GetScalars("Q");

    tetrahedron_vector.reserve(number_of_cells);

    for (size_t k = 0; k < number_of_cells; k++) {
        /*
         * in newer vtk versions should be replaced with range based iterators.
         * which are already in nightly releases 04.11.19
         */
        vtkSmartPointer<vtkPoints> points = unstructured_grid->GetCell(k)->GetPoints();
        std::array<std::array<double, 3>, 4> tmp_points{};
        for (size_t i = 0; i < 4; i++) {
            double *p = points->GetPoint(i);
            std::copy(p, p + 3, tmp_points[i].begin());
        }

        double *tmp_q = scalars_q->GetTuple(k);
        double *tmp_alpha = scalars_alpha->GetTuple(k);
        tetra im(tmp_points, *tmp_alpha, *tmp_q);

        tetrahedron_vector.push_back(im);
    }
    return tetrahedron_vector;
}

std::vector<tetra> app::get_tetrahedron_vector_binary(const std::string& filename) {
    std::vector<tetra> tetrahedron_vector{};
    vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid{init_vtk_grid(filename)};
    vtkSmartPointer<vtkCellData> scalar_data = unstructured_grid->GetCellData();

    size_t number_of_cells = unstructured_grid->GetNumberOfCells();

//    scalar_data->Print(std::cout);
    vtkSmartPointer<vtkDataArray> scalars_alpha = scalar_data->GetScalars("AbsorpCoef");
    vtkSmartPointer<vtkDataArray> scalars_q = scalar_data->GetScalars("radEnLooseRate");

    tetrahedron_vector.reserve(number_of_cells);

    for (size_t k = 0; k < number_of_cells; k++) {
        /*
         * in newer vtk versions should be replaced with range based iterators.
         * which are already in nightly releases 04.11.19
         */
        vtkSmartPointer<vtkPoints> points = unstructured_grid->GetCell(k)->GetPoints();
        std::array<std::array<double, 3>, 4> tmp_points{};
        for (size_t i = 0; i < 4; i++) {
            double *p = points->GetPoint(i);
            std::copy(p, p + 3, tmp_points[i].begin());
        }

        double *tmp_q = scalars_q->GetTuple(k);
        double *tmp_alpha = scalars_alpha->GetTuple(k);
        tetra im(tmp_points, *tmp_alpha, *tmp_q);

        tetrahedron_vector.push_back(im);
    }
    return tetrahedron_vector;
}

void app::rotate_tetrahedron_vector(std::vector<tetra>& tetrahedron_vector, double angle) {
    const size_t tetra_vec_size = tetrahedron_vector.size();

#pragma omp parallel for default(none) shared(tetrahedron_vector, angle) schedule(dynamic, 8)
    for (size_t i = 0; i < tetra_vec_size; i++) {
        tetrahedron_vector[i].rotate_x(angle);
    }
}

/*
 * global domain boundaries
 * {x_max, x_min, y_max, y_min}
 */
std::array<double, 4> app::get_domain_boundaries(std::vector<tetra>& tetrahedron_vector) {
    auto pre_tmp = tetrahedron_vector.at(0).get_boundaries();

    double x_max = pre_tmp[0];
    double x_min = pre_tmp[1];
    double y_max = pre_tmp[2];
    double y_min = pre_tmp[3];

    const size_t size_of_vec = tetrahedron_vector.size();

    /*
     * TODO:
     * 1. move to main!
     */
    omp_set_num_threads(AMOUNT_OF_THREADS);
#pragma omp parallel for default(none) shared(tetrahedron_vector) reduction(max: x_max, y_max) reduction(min: x_min, y_min) schedule(dynamic, 8)
    for (size_t i = 0; i < size_of_vec; i++) {
        auto tmp = tetrahedron_vector[i].get_boundaries();
        if (tmp[0] > x_max) {
            x_max = tmp[0];
        }
        if (tmp[1] < x_min) {
            x_min = tmp[1];
        }
        if (tmp[2] > y_max) {
            y_max = tmp[2];
        }
        if (tmp[3] < y_min) {
            y_min = tmp[3];
        }
    }

    std::array<double, 4> res{x_max, x_min, y_max, y_min};
//    std::cout << x_max << " " << x_min << " " << y_max << " " << y_min << std::endl;

    return res;
}

void app::print_boundaries(const std::array<double, 4>& global_boundaries) {
    std::cout
        << "plane boundaries" << std::endl
        << "x_max = " << global_boundaries[0] << std::endl
        << "x_min = " << global_boundaries[1] << std::endl
        << "y_max = " << global_boundaries[2] << std::endl
        << "y_min = " << global_boundaries[3] << std::endl;
}

plane app::init_plane_grid(size_t res_x, size_t res_y, std::array<double, 4>& global_boundaries) {
    return plane(res_x, res_y, global_boundaries);
}

void app::find_tetrahedron_vector_intersections_with_lines(const std::vector<tetra>& tetrahedron_vector,
                                                           plane& task_plane) {
    auto t1 = std::chrono::high_resolution_clock::now();
    task_plane.find_intersections_with_tetra_vector(tetrahedron_vector);
    auto t2 = std::chrono::high_resolution_clock::now();

    auto timer = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    std::cout << "find all intersections complete in: " << timer << " ms" << std::endl;
}

std::pair<float_matrix, float_matrix>
app::trace_rays(plane& current_plane, std::vector<tetra>& tetrahedron_vector, tetra_value value_alpha,
                tetra_value value_q) {
    auto t1 = std::chrono::high_resolution_clock::now();
//    std::cout << current_plane.count_all_intersections() << std::endl;

    auto res = current_plane.trace_rays(tetrahedron_vector, value_alpha, value_q);

    auto t2 = std::chrono::high_resolution_clock::now();
    auto timer = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();

    std::cout << "trace all rays complete in: " << timer << " ms" << std::endl;

    return res;
}

void
app::produce_result(std::pair<float_matrix, float_matrix>& data, size_t res_x, size_t res_y,
                    const std::string& result_filename) {
    vtkSmartPointer<vtkImageData> imageData =
        vtkSmartPointer<vtkImageData>::New();
    imageData->SetDimensions(res_x, res_y, 1);
    imageData->AllocateScalars(VTK_DOUBLE, 2);

    int *dims = imageData->GetDimensions();

    for (int y = 0; y < dims[1]; y++) {
        for (int x = 0; x < dims[0]; x++) {
            auto *pixel = static_cast<double *>(imageData->GetScalarPointer(x, y, 0));
            pixel[0] = data.first[x][y];
            pixel[1] = data.second[x][y];
        }
    }

    vtkSmartPointer<vtkXMLImageDataWriter> writer =
        vtkSmartPointer<vtkXMLImageDataWriter>::New();
    writer->SetFileName(result_filename.c_str());
    writer->SetInputData(imageData);
    writer->Write();
}
