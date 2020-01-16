#pragma once

#include <object3d_base.hpp>
#include <string>

class object3d_accretion_disk : public object3d_base {
public:
    explicit object3d_accretion_disk(const std::string& filename);
};