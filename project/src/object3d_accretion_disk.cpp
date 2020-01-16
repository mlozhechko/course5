#include <object3d_accretion_disk.hpp>

object3d_accretion_disk::object3d_accretion_disk(const std::string& filename) {
    object3d_base::read_vtk_file(filename, {"AbsorpCoef", "radEnLooseRate"});
}
