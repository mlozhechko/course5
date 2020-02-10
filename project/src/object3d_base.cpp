#include <object3d_base.hpp>

static vtkSmartPointer<vtkUnstructuredGrid> init_vtk_grid(const std::string& filename) {
    auto reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader->SetFileName(filename.c_str());
    reader->SetReadAllScalars(true);
    reader->Update();

    return reader->GetOutput();
}

void object3d_base::read_vtk_file(const std::string& filename,
                                  const std::vector<std::string>& scalar_labels) {
    vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid{init_vtk_grid(filename)};
    vtkSmartPointer<vtkCellData> scalar_data = unstructured_grid->GetCellData();

    size_t number_of_cells = unstructured_grid->GetNumberOfCells();

    std::array<vtkSmartPointer<vtkDataArray>, 2> scalars{};
    const size_t amount_of_scalars = scalar_labels.size();

    for (size_t i = 0; i < amount_of_scalars; i++) {
        scalars[i] = scalar_data->GetScalars(scalar_labels[i].data());
    }

    _data->reserve(number_of_cells);

    for (size_t k = 0; k < number_of_cells; k++) {
        /*
         * in newer vtk versions should be replaced with range based iterators.
         * which are already in nightly releases 04.11.19
         */
        vtkSmartPointer<vtkPoints> points = unstructured_grid->GetCell(k)->GetPoints();
        std::array<std::array<double, 3>, 4> tmp_points{};
        for (size_t i = 0; i < 4; i++) {
            double* p = points->GetPoint(i);
            std::copy(p, p + 3, tmp_points[i].begin());
        }

        std::array<double, 2> scalar_values{};
        for (size_t i = 0; i < amount_of_scalars; i++) {
            scalar_values[i] = *(scalars[i]->GetTuple(k));
        }

        tetra im(tmp_points, scalar_values[0], scalar_values[1], tetra_type::transparent);
        _data->push_back(im);
    }
}

static std::array<double, 3> rotate_vector_around_y_axis(const std::array<double, 3>& line_vector,
                                                         double alpha) {
    std::array<double, 3> res{};
    res[0] = line_vector[0] * cos(alpha) + line_vector[2] * sin(alpha);
    res[1] = line_vector[1];
    res[2] = -line_vector[0] * sin(alpha) + line_vector[2] * cos(alpha);

    return res;
}

static std::array<double, 3> rotate_vector_around_z_axis(const std::array<double, 3>& line_vector,
                                                         double alpha) {
    std::array<double, 3> res{};
    res[0] = line_vector[0] * cos(alpha) + line_vector[1] * sin(alpha);
    res[1] = -line_vector[0] * sin(alpha) + line_vector[1] * cos(alpha);
    res[2] = line_vector[2];

    return res;
}

static void add_vector(std::array<double, 3>& arr1, const std::array<double, 3>& arr2) {
    for (size_t i = 0; i < arr1.size(); i++) {
        arr1[i] += arr2[i];
    }
}

void object3d_base::init_polar(
    const std::function<double(std::array<double, 3>)>& potential_function, double x0, double y0,
    double z0, double level_value, double step, double angle_step, tetra_type arg_tetra_type,
    double tetra_v1, double tetra_v2) {
    const std::array<double, 3> def_step_vector{0.001, 0, 0};

    double step_angle_x = PI / angle_step;
    double step_angle_y = PI / angle_step;

    double angle_x = 0;
    double angle_y = -PI + step_angle_y;

    std::array<double, 3> bottom_point{};
    std::array<double, 3> top_point{};

    {
        std::array<double, 3> point_trace_vec{x0, y0, z0};
        std::array<double, 3> step_vector{0, 0, step};

        double f_value{};
        do {
            add_vector(point_trace_vec, step_vector);
            //            f_value = roche_lobe_potential(point_trace_vec, ACC_X0, donor_pos_x,
            //            mass_center_pos_x);
            f_value = potential_function(point_trace_vec);
        } while (f_value < level_value);

        top_point = point_trace_vec;
        point_trace_vec = {x0, y0, z0};
        step_vector = {0, 0, -step};

        do {
            add_vector(point_trace_vec, step_vector);
            f_value = potential_function(point_trace_vec);
        } while (f_value < level_value);

        bottom_point = point_trace_vec;
    }

    std::vector<std::vector<std::array<double, 3>>> polar_points{};

    while (angle_y < (PI - step_angle_y + std::numeric_limits<double>::epsilon())) {
        std::array<double, 3> y_step_vector = rotate_vector_around_z_axis(def_step_vector, angle_y);
        angle_y += step_angle_y;

        std::vector<std::array<double, 3>> line{};
        while (angle_x < 2 * PI - step_angle_x + std::numeric_limits<double>::epsilon()) {
            std::array<double, 3> point_trace_vec{x0, y0, z0};
            std::array<double, 3> step_vector = rotate_vector_around_y_axis(y_step_vector, angle_x);

            double f_value{};
            do {
                add_vector(point_trace_vec, step_vector);
                f_value = potential_function(point_trace_vec);
            } while (f_value < level_value);

            line.push_back(point_trace_vec);
            angle_x += step_angle_x;
        }

        polar_points.push_back(std::move(line));
        angle_x = 0;
    }

    /*
     * create vector of tetras of lobe
     *
     * main cycle from bottom to top
     */

    const std::array<double, 3> center{x0, y0, z0};
    std::vector<tetra> object3d_data{};

    const size_t line_len = polar_points[0].size();
    for (size_t i = 1; i < line_len; i++) {
        std::array<std::array<double, 3>, 4> points_bot{center, bottom_point, polar_points[0][i],
                                                        polar_points[0][i - 1]};
        object3d_data.emplace_back(points_bot, tetra_v1, tetra_v2, arg_tetra_type);
    }
    std::array<std::array<double, 3>, 4> points_bot{center, bottom_point, polar_points[0][0],
                                                    polar_points[0][line_len - 1]};
    object3d_data.emplace_back(points_bot, tetra_v1, tetra_v2, arg_tetra_type);

    const size_t h_len = polar_points.size() - 1;
    for (size_t i = 1; i < line_len; i++) {
        std::array<std::array<double, 3>, 4> points_top{center, top_point, polar_points[h_len][i],
                                                        polar_points[h_len][i - 1]};
        object3d_data.emplace_back(points_top, tetra_v1, tetra_v2, arg_tetra_type);
    }

    std::array<std::array<double, 3>, 4> points_top{center, top_point, polar_points[h_len][0],
                                                    polar_points[0][line_len - 1]};
    object3d_data.emplace_back(points_top, tetra_v1, tetra_v2, arg_tetra_type);

    for (size_t i = 1; i < (h_len + 1); i++) {
        for (size_t j = 1; j < line_len; j++) {
            std::array<std::array<double, 3>, 4> points1{
                center, polar_points[i - 1][j - 1], polar_points[i - 1][j], polar_points[i][j - 1]};
            object3d_data.emplace_back(points1, tetra_v1, tetra_v2, arg_tetra_type);

            std::array<std::array<double, 3>, 4> points2{
                center, polar_points[i][j - 1], polar_points[i][j], polar_points[i - 1][j]};
            object3d_data.emplace_back(points2, tetra_v1, tetra_v2, arg_tetra_type);
        }

        std::array<std::array<double, 3>, 4> points1{center, polar_points[i - 1][line_len - 1],
                                                     polar_points[i - 1][0],
                                                     polar_points[i][line_len - 1]};
        object3d_data.emplace_back(points1, tetra_v1, tetra_v2, arg_tetra_type);

        std::array<std::array<double, 3>, 4> points2{center, polar_points[i][line_len - 1],
                                                     polar_points[i][0], polar_points[i - 1][0]};
        object3d_data.emplace_back(points2, tetra_v1, tetra_v2, arg_tetra_type);
    }

    *_data = std::move(object3d_data);
}

std::shared_ptr<std::vector<tetra>> object3d_base::get_pointer() {
    return _data;
}

void object3d_base::rotate_around_x_axis(double angle) {
    const size_t tetra_vec_size = _data->size();

#pragma omp parallel for default(none) shared(_data, angle) schedule(dynamic, 8)
    for (size_t i = 0; i < tetra_vec_size; i++) {
        (*_data)[i].rotate_around_x_axis(angle);
    }
}

void object3d_base::rotate_around_y_axis(double angle, double x0) {
    const size_t tetra_vec_size = _data->size();

#pragma omp parallel for default(none) shared(_data, angle, x0) schedule(dynamic, 8)
    for (size_t i = 0; i < tetra_vec_size; i++) {
        (*_data)[i].rotate_around_y_axis(angle, x0);
    }
}

std::array<double, 4> object3d_base::get_boundaries() {
    auto pre_tmp = _data->at(0).get_boundaries();

    double x_max = pre_tmp[0];
    double x_min = pre_tmp[1];
    double y_max = pre_tmp[2];
    double y_min = pre_tmp[3];

    const size_t size_of_vec = _data->size();

    /*
     * TODO:
     * 1. move to main!
     */
    //        omp_set_num_threads(AMOUNT_OF_THREADS);

#pragma omp parallel for default(none) shared(_data) reduction(max                                 \
                                                               : x_max, y_max)                     \
    reduction(min                                                                                  \
              : x_min, y_min) schedule(dynamic, 8)
    for (size_t i = 0; i < size_of_vec; i++) {
        auto tmp = (*_data)[i].get_boundaries();
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

    return std::array<double, 4>{x_max, x_min, y_max, y_min};
}