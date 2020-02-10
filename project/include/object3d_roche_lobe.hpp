#pragma once

#include <object3d_base.hpp>

class object3d_roche_lobe : public object3d_base {
public:
    object3d_roche_lobe(const point& pos_accretor, double L, double donor_angle_around_y,
                        double m_accretor, double m_donor, double def_omega);
};