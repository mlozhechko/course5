#include <iostream>
#include <string>

#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>

constexpr double pi = 3.14159265358979323846;

class lite_tetrahedron {
public:
    explicit lite_tetrahedron(const std::array<std::array<double, 3>, 4>& points, double alpha, double q)
    : _q(q), _alpha(alpha), _points(points){}

    const std::array<double, 3>& operator[](size_t i) const {
        return _points[i];
    }

    /*
     * rotation by angle by y axis
     */
    void rotate_x(double angle) {
        for (auto& i: _points) {
            point_rotate_x(angle, i);
        }
    }

    /*
     * get tetrahedrons boundaries
     * std::array<double, 4> = {x_max, x_min, y_max, y_min}
     */
    std::array<double, 4> get_boundaries() {
        /*
         * TODO:
         * 1. rewrite with points instead of double
         */
        std::array<double, 4> res{};
        res[1] = res[0] = _points[0][0];
        res[3] = res[2] = _points[0][1];

        for (size_t i = 1; i < 4; i++) {
            /*
             * GPU optimization to avoid if statements (blank)
             *
             * TODO:
             * check for accuracy losses
             */
            auto flag = static_cast<double>(res[0] < _points[i][0]);
            res[0] = flag * _points[i][0] + (1 - flag) * res[0];

            flag = static_cast<double>(res[1] >_points[i][0]);
            res[1] = flag * _points[i][0] + (1 - flag) * res[1];

            flag = static_cast<double>(res[2] < _points[i][1]);
            res[2] = flag * _points[i][1] + (1 - flag) * res[2];

            flag = static_cast<double>(res[3] > _points[i][1]);
            res[3] = flag * _points[i][1] + (1 - flag) * res[3];
        }

        return res;
    }

private:
    static void point_rotate_x(double alpha, std::array<double, 3>& point) {
        double tmp_point1 = point[1];
        point[1] = point[1] * cos(alpha) - point[2] * sin(alpha);
        point[2] = tmp_point1 * sin(alpha) + point[2] * cos(alpha);
    }

    friend std::ostream& operator<<(std::ostream& os, const lite_tetrahedron& lt);

    std::array<std::array<double, 3>, 4> _points{};
    double _q;
    double _alpha;
};

std::ostream& operator<<(std::ostream& os, const lite_tetrahedron& lt)
{
    os << "lite tetrahedron. q = " << lt._q << ", alpha = " << lt._alpha << std::endl;
    for (const auto& it: lt._points) {
        os << it[0] << " " << it[1] << " " << it[2] << std::endl;
    }
    return os;
}

class line {
public:
    explicit line(double x, double y) : _x(x), _y(y) {};

    double x() {
        return _x;
    }
    double y() {
        return _y;
    }
    void add_possible_tetra_id(size_t id) {
        _tetra_id.push_back(id);
    }
private:
    double _x, _y;
    std::vector<size_t> _tetra_id;
};

class plane {
public:
    plane(size_t res_x, size_t res_y, const std::array<double, 4>& global_boundaries) {
        double delta_x(global_boundaries[0] - global_boundaries[1]);
        double delta_y(global_boundaries[2] - global_boundaries[3]);

        _global_boundaries = global_boundaries;

        _step_x = delta_x / (res_x - 1.);
        _step_y = delta_y / (res_y - 1.);

        _lines.resize(res_x);
        double current_x = global_boundaries[1]; //x_min
        for (size_t i = 0; i < res_x; i++) {
            double current_y = global_boundaries[3]; //y_min
            _lines.at(i).reserve(res_y);
            for (size_t j = 0; j < res_y; j++) {
                _lines.at(i).push_back(line(current_x, current_y));
                current_y = current_y + _step_y;
            }
            current_x = current_x + _step_x;
        }
    }

    /*
     * return amount of lines assumed to intersect tetrahedron
     */
    size_t add_possible_tetra_from_boundaries(const std::array<double, 4>& boundaries, size_t id) {
        /*
         * TODO:
         * 1. rewrite method to add all explicit intersection tetras from boundaries
         */

        /*
         * x_start = (x_min - global_x_min) / step
         * */
        size_t x_start = ceil((boundaries[1] - _global_boundaries[1]) / _step_x);

        /*
         * x_end = (x_max - global_x_min) / step
         * */
        size_t x_end = ceil((boundaries[0] - _global_boundaries[1]) / _step_x);

        /*
         * y_start = (y_min - global_y_min) / step
         */
        size_t y_start = ceil((boundaries[3] - _global_boundaries[3]) / _step_y);

        /*
         * y_start = (y_min - global_y_min) / step
         */
        size_t y_end = ceil((boundaries[2] - _global_boundaries[3]) / _step_y);

        for (size_t i = x_start; i < x_end; i++) {
            for (size_t j = y_start; j < y_end; j++) {
                _lines.at(i).at(j).add_possible_tetra_id(id);
            }
        }

        return (x_end - x_start) * (y_end - y_start);
    }

private:
    /*
     * [] sequence is [x][y]
     */
    std::vector<std::vector<line>> _lines;
    std::array<double, 4> _global_boundaries;

    double _step_x, _step_y;

};

vtkSmartPointer<vtkUnstructuredGrid> init_grid(const std::string& filename) {
    auto reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader->SetFileName(filename.c_str());
    reader->SetReadAllScalars(true);
    reader->Update();

//    reader->Print(std::cout);
    return reader->GetOutput();
}

std::vector<lite_tetrahedron> init_tetrahedron_vector(/* some boost program_options arguments */) {
    std::vector<lite_tetrahedron> tetrahedron_vector{};
    vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid{init_grid("test.vtk")};
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

int main() {
    std::vector<lite_tetrahedron> tetrahedron_vector = std::move(init_tetrahedron_vector());
    std::cout << tetrahedron_vector[24000] << std::endl;

    double angle = pi / 12.;
    for (auto& i: tetrahedron_vector) {
        i.rotate_x(angle);
    }

    std::cout << tetrahedron_vector[1000] << std::endl;

    std::vector<std::array<double, 4>> tetrahedron_boundaries{};
    tetrahedron_boundaries.reserve(tetrahedron_vector.size());

    /*
     * global domain boundaries
     * {x_max, x_min, y_max, y_min}
     */
    std::array<double, 4> global_boundaries{tetrahedron_vector.at(0).get_boundaries()};
    std::transform(tetrahedron_vector.begin(), tetrahedron_vector.end(), std::back_inserter(tetrahedron_boundaries), [&global_boundaries](auto& x) {
        /*
         * TODO:
         * check for bottleneck
         */
        auto tmp = x.get_boundaries();

        /*
         * init global boundaries
         */
        global_boundaries[0] = std::max(global_boundaries[0], tmp[0]);
        global_boundaries[1] = std::min(global_boundaries[1], tmp[1]);

        global_boundaries[2] = std::max(global_boundaries[2], tmp[2]);
        global_boundaries[3] = std::min(global_boundaries[3], tmp[3]);
        return tmp;
    });

    std::cout << tetrahedron_boundaries[1000][3] << std::endl;

    std::cout
        << "plane boundaries"
        << "x_max = " << global_boundaries[0] << std::endl
        << "x_min = " << global_boundaries[1] << std::endl
        << "y_max = " << global_boundaries[2] << std::endl
        << "y_min = " << global_boundaries[3] << std::endl;

    /*
     * TODO
     * create program arguments for resolution
     */
    size_t grid_resolution_x(620), grid_resolution_y(480);
    plane current_plane(grid_resolution_x, grid_resolution_y, global_boundaries);

    size_t min_intersections(0), max_intersections(0);
    min_intersections = max_intersections = current_plane.add_possible_tetra_from_boundaries(tetrahedron_boundaries[0], 0);

    for (size_t i = 1; i < tetrahedron_boundaries.size(); i++) {
        size_t res = current_plane.add_possible_tetra_from_boundaries(tetrahedron_boundaries[i], i);
        min_intersections = std::min(res, min_intersections);
        max_intersections = std::max(res, max_intersections);
    }

    std::cout
        << "max possible intersections = " << max_intersections << std::endl
        << "min possible intersections = " << min_intersections << std::endl;

    return 0;
}