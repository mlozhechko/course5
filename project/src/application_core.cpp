#include <application_core.hpp>

application_core& application_core::instance() {
    static application_core instance;
    return instance;
}

application_core::application_core() = default;

vtkSmartPointer<vtkUnstructuredGrid> application_core::init_vtk_grid(const std::string& tmp_filename) {
    auto reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader->SetFileName(tmp_filename.c_str());
    reader->SetReadAllScalars(true);
    reader->Update();

//    reader->Print(std::cout);
    return reader->GetOutput();
}

std::vector<lite_tetrahedron> application_core::get_tetrahedron_vector(/* some boost program_options arguments */) {
    std::vector<lite_tetrahedron> tetrahedron_vector{};
    vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid{init_vtk_grid(_filename)};
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

void application_core::rotate_tetrahedron_vector(std::vector <lite_tetrahedron>& tetrahedron_vector, double angle) {
    for (auto& i: tetrahedron_vector) {
        i.rotate_x(angle);
    }
}

/*
 * global domain boundaries
 * {x_max, x_min, y_max, y_min}
 */
std::array<double, 4> application_core::get_global_boundaries(std::vector<lite_tetrahedron>& tetrahedron_vector) {
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

void application_core::print_boundaries(const std::array<double, 4>& global_boundaries) {
    std::cout
        << "plane boundaries" << std::endl
        << "x_max = " << global_boundaries[0] << std::endl
        << "x_min = " << global_boundaries[1] << std::endl
        << "y_max = " << global_boundaries[2] << std::endl
        << "y_min = " << global_boundaries[3] << std::endl;
}

plane application_core::init_plane_grid(std::array<double, 4>& global_boundaries) {
    size_t grid_resolution_x(application_core::_res_x), grid_resolution_y(application_core::_res_y);
    return plane(grid_resolution_x, grid_resolution_y, global_boundaries);
}

void application_core::set_filename(const std::string& filename) {
    _filename = filename;
}

void application_core::set_resolution(const std::array<size_t, 2>& res_arr) {
    _res_x = res_arr[0];
    _res_y = res_arr[1];
}