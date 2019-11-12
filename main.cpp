#include <iostream>
#include <string>
#include <bitset>
#include <algorithm>

#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <chrono>

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
         * 1. remove(!)
         */
        std::array<double, 4> res{};
        res[1] = res[0] = _points[0][0];
        res[3] = res[2] = _points[0][1];

        for (size_t i = 1; i < 4; i++) {
            /*
             * GPU optimization to avoid if statements (blank)
             *
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
    explicit line(double x, double y) : _x(x), _y(y) {
        _tetra_intersections.reserve(5000);
    };
    double x() {
        return _x;
    }
    double y() {
        return _y;
    }
    void add_tetra_intersection(size_t id, const std::bitset<4>& polygon_intersections) {
//        std::cout << "{" << x() << "," << y() << "}," << polygon_intersections << std::endl;
        _tetra_intersections.emplace_back(id, polygon_intersections);
    }
private:
    double _x, _y;
    std::vector<std::pair<size_t, std::bitset<4>>> _tetra_intersections;
};

class plane {
public:
    plane(size_t res_x, size_t res_y, const std::array<double, 4>& global_boundaries)
    : _global_boundaries(global_boundaries) {
        double delta_x(global_boundaries[0] - global_boundaries[1]);
        double delta_y(global_boundaries[2] - global_boundaries[3]);

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
     * return amount of lines to intersect polygons
     */
    size_t find_intersections_with_tetrahedron(const lite_tetrahedron& tetra, size_t id) {
        /*
         * 0001 0, 1, 2 points
         * 0010 0, 1, 3 points
         * 0100 0, 2, 3 points
         * 1000 1, 2, 3 points
         */

        size_t sum(0);
        sum += find_intersections_with_polygon({tetra[0], tetra[1], tetra[2]}, id, 0b1110);
        sum += find_intersections_with_polygon({tetra[0], tetra[1], tetra[3]}, id, 0b1101);
        sum += find_intersections_with_polygon({tetra[0], tetra[2], tetra[3]}, id, 0b1011);
        sum += find_intersections_with_polygon({tetra[1], tetra[2], tetra[3]}, id, 0b0111);

        if (sum % 2 == 1) {
            throw std::runtime_error("critical error. odd number of intersections");
        }

        return sum / 2;
    }

private:
    /*
     * all line equations is on X,Y plane
     */

    static inline double line_common_eq(std::array<double, 3> p1, std::array<double, 3> p2, std::array<double, 3> pos) {
         return (p2[1] - p1[1]) * pos[0] + (p1[0] - p2[0]) * pos[1] + (p2[0] * p1[1] - p1[0] * p2[1]);
    }

    static inline double line_rev_function_eq(std::array<double, 3> p1, std::array<double, 3> p2, double y) {
        return (p1[0] - p2[0]) * (y - p1[1]) / (p1[1] - p2[1]) + p1[0];
    }

    size_t find_intersections_with_polygon(std::array<std::array<double, 3>, 3> points, size_t id, std::bitset<4> polygon) {
        size_t counter(0);

        std::sort(points.begin(), points.end(), [](auto a, auto b){
            return a[1] > b[1];
        });

        /*
         * find p[1] pos relative to p[0] to p[2] line
         */
        double rel_qual = line_common_eq(points[0], points[2], points[1]);

        /*
         * let L be line passing through points p[0] and p[2]
         *
         * check if L is ascending and p[1] is below the line
         */
        bool asc_above = (points[0][0] >= points[2][0]) && (rel_qual >= 0);

        /*
         * check if L is descending and p[1] is above the line
         */
        bool des_below = (points[0][0] < points[2][0]) && (rel_qual > 0);

        /*
         * p[1] may on the left side or on right side of the line
         * also p[1].y is greater than p[2].y and lower than p[0].y
         *
         * therefore we can fully describe position of triangle of plane, and use it
         * to find grid area that intersects polygon.
         *
         * position is true if p[2] is on the right side of L
         */
        bool position = !(asc_above || des_below);

        double y_max = points[0][1];
        double y_min = points[2][1];

        size_t y_max_index = std::floor((y_max - _global_boundaries[3]) / _step_y);
        size_t y_min_index = std::ceil((y_min - _global_boundaries[3]) / _step_y);

        size_t y_index_it = y_min_index;
        double y_it = _lines[0][y_index_it].y();

        /*
         * iterate over triangle's edges and find intersection points
         */
        for (; y_index_it <= y_max_index; y_index_it++) {
            double x_max(0), x_min(0);
            if (position) {
                x_min = line_rev_function_eq(points[0], points[2], y_it);
                if (y_it < points[1][1]) {
                    x_max = line_rev_function_eq(points[2], points[1], y_it);
                } else {
                    x_max = line_rev_function_eq(points[0], points[1], y_it);
                }
            } else {
                x_max = line_rev_function_eq(points[0], points[2], y_it);
                if (y_it < points[1][1]) {
                    x_min = line_rev_function_eq(points[2], points[1], y_it);
                } else {
                    x_min = line_rev_function_eq(points[0], points[1], y_it);
                }
            }

            size_t x_max_index = std::floor((x_max - _global_boundaries[1]) / _step_x);
            size_t x_min_index = std::ceil((x_min - _global_boundaries[1]) / _step_x);

            for (size_t i = x_min_index; i <= x_max_index; i++) {
                _lines[i][y_index_it].add_tetra_intersection(id, polygon);
                counter++;
            }

            y_it = y_it + _step_y;
        }

        return counter;
    }

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

    /*
     * global domain boundaries
     * {x_max, x_min, y_max, y_min}
     */
    std::array<double, 4> global_boundaries{tetrahedron_vector.at(0).get_boundaries()};
    std::for_each(tetrahedron_vector.begin(), tetrahedron_vector.end(), [&global_boundaries](auto& x) {
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
    });

    std::cout
        << "plane boundaries" << std::endl
        << "x_max = " << global_boundaries[0] << std::endl
        << "x_min = " << global_boundaries[1] << std::endl
        << "y_max = " << global_boundaries[2] << std::endl
        << "y_min = " << global_boundaries[3] << std::endl;

    /*
     * TODO
     * create program arguments for resolution
     */
    size_t grid_resolution_x(1920), grid_resolution_y(1080);
    plane current_plane(grid_resolution_x, grid_resolution_y, global_boundaries);

    size_t min_intersections(0), max_intersections(0);
    min_intersections = max_intersections = current_plane.find_intersections_with_tetrahedron(tetrahedron_vector[0], 0);


    auto t1 = std::chrono::high_resolution_clock::now();

    for (size_t i = 1; i < tetrahedron_vector.size(); i++) {
        size_t res = current_plane.find_intersections_with_tetrahedron(tetrahedron_vector[i], i);
        min_intersections = std::min(res, min_intersections);
        max_intersections = std::max(res, max_intersections);
    }

    auto t2 = std::chrono::high_resolution_clock::now();
    auto timer = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    std::cout << "find all intersections complete in: " << timer << " ms" << std::endl;

    std::cout << max_intersections << std::endl;
    std::cout << min_intersections << std::endl;

    return 0;
}