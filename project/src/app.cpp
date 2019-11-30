#include <app.hpp>


vtkSmartPointer<vtkUnstructuredGrid> app::init_vtk_grid(const std::string& filename) {
    auto reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader->SetFileName(filename.c_str());
    reader->SetReadAllScalars(true);
    reader->Update();

//    reader->Print(std::cout);
    return reader->GetOutput();
}

std::vector<lite_tetrahedron> app::get_tetrahedron_vector(const std::string& filename) {
    std::vector<lite_tetrahedron> tetrahedron_vector{};
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
        lite_tetrahedron im(tmp_points, *tmp_alpha, *tmp_q);

        tetrahedron_vector.push_back(im);
    }
    return tetrahedron_vector;
}

void app::rotate_tetrahedron_vector(std::vector <lite_tetrahedron>& tetrahedron_vector, double angle) {
    for (auto& i: tetrahedron_vector) {
        i.rotate_x(angle);
    }
}

/*
 * global domain boundaries
 * {x_max, x_min, y_max, y_min}
 */
std::array<double, 4> app::get_domain_boundaries(std::vector<lite_tetrahedron>& tetrahedron_vector) {
    std::array<double, 4> global_boundaries{tetrahedron_vector.at(0).get_boundaries()};
    std::for_each(tetrahedron_vector.begin(), tetrahedron_vector.end(), [&global_boundaries](auto& x) {
        auto tmp = x.get_boundaries();
        global_boundaries[0] = std::max(global_boundaries[0], tmp[0]);
        global_boundaries[1] = std::min(global_boundaries[1], tmp[1]);
        global_boundaries[2] = std::max(global_boundaries[2], tmp[2]);
        global_boundaries[3] = std::min(global_boundaries[3], tmp[3]);
    });

    return global_boundaries;
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

void app::find_tetrahedron_vector_intersections_with_lines(const std::vector<lite_tetrahedron>& tetrahedron_vector,
                                                           plane& task_plane) {
    size_t min_intersections{0}, max_intersections{0};
    min_intersections = max_intersections = task_plane.find_intersections_with_tetrahedron(tetrahedron_vector[0], 0);

    auto t1 = std::chrono::high_resolution_clock::now();

    for (size_t i = 1; i < tetrahedron_vector.size(); i++) {
        size_t res = task_plane.find_intersections_with_tetrahedron(tetrahedron_vector[i], i);
        min_intersections = std::min(res, min_intersections);
        max_intersections = std::max(res, max_intersections);
    }

    auto t2 = std::chrono::high_resolution_clock::now();

    auto timer = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    std::cout
        << "find all intersections complete in: " << timer << " ms"
//        << std::endl
//        << "maximum and minimum amound of_intersections with tetrahedron "
//        << max_intersections << " " << min_intersections
        << std::endl;
}

std::pair<float_matrix, float_matrix>
app::trace_rays(plane& current_plane, std::vector<lite_tetrahedron>& tetrahedron_vector, tetra_value value_alpha,
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

    int* dims = imageData->GetDimensions();

    for (int y = 0; y < dims[1]; y++)
    {
        for (int x = 0; x < dims[0]; x++)
        {
            auto* pixel = static_cast<double*>(imageData->GetScalarPointer(x, y,0));
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
