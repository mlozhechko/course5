#include <object3d_roche_lobe.hpp>

static double vector_2_norm(const std::array<double, 3>& vec) {
    double sum{};
    for (const auto& i : vec) {
        sum += i * i;
    }
    return sqrt(sum);
}

static std::array<double, 3> vector_multiplication(const std::array<double, 3>& m1, const std::array<double, 3>& m2) {
    std::array<double, 3> res{};
    res[0] = m1[1] * m2[2] - m1[2] * m2[1];
    res[1] = -(m1[0] * m2[2]) + (m1[2] * m2[0]);
    res[2] = m1[0] * m2[1] - m1[1] * m2[0];

    return res;
}

object3d_roche_lobe::object3d_roche_lobe(const point& pos_accretor, double dist, double donor_angle_around_y,
                                         double m_accretor, double m_donor, double def_omega) {
    const double donor_pos_x = pos_accretor[0] - dist;
    const double mass_center_pos_x = (donor_pos_x * m_donor + pos_accretor[0] * m_accretor) / (m_accretor + m_donor);
    const double alpha = m_donor / (m_accretor + m_donor);
    /* const */ double lagrange1_pos_x = mass_center_pos_x - dist * (1. - pow((alpha / 3.), 1. / 3.));

    /*
     * for some reason defined in article lagrange_pos_x is another value
     */
    lagrange1_pos_x = 0.35515;

    auto roche_lobe_potential_func = [&](const std::array<double, 3>& r) -> double {
        double acc_denominator = vector_2_norm({r[0] - pos_accretor[0], r[1], r[2]});
        //    std::cout << "acc denominator: " << acc_denominator << std::endl;
        double donor_denominator = vector_2_norm({r[0] - donor_pos_x, r[1], r[2]});
        //    std::cout << "donor denominator: " << donor_denominator << std::endl;
        double tmp_omega_mult =
            vector_2_norm(vector_multiplication({r[0] - mass_center_pos_x, r[1], r[2]}, {0, def_omega, 0}));
        double omega = (1. / 2.) * tmp_omega_mult * tmp_omega_mult;
        //    std::cout << "omega: " << omega << std::endl;

        double F = -((G_SOL * m_accretor) / acc_denominator) - ((G_SOL * m_donor) / donor_denominator) - omega;
        return F;
    };

    const double lagrange_potential = roche_lobe_potential_func({lagrange1_pos_x, 0, 0});
    init_polar(roche_lobe_potential_func, donor_pos_x, 0, 0, lagrange_potential, 0.001, 128, tetra_type::solid, app::instance().config.roche_lobe_solid_color);
    rotate_around_y_axis(donor_angle_around_y, ACC_X0);
}
