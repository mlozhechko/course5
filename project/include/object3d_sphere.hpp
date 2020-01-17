#pragma once

#include <object3d_base.hpp>

class object3d_sphere : public object3d_base {
public:
    object3d_sphere(const point& center, double R);
};