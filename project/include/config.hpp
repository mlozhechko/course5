#pragma once

/*
 * Render configuration
 */

class app{
public:

    static app& instance() {
        static app instance;
        return instance;
    };

    struct config_str {
        std::string file = {};
        std::string destination = {};
        int threads = 1;
        size_t resolution_x = 1100;
        size_t resolution_y = 900;
        double angle_around_x = 0;
        double angle_around_y = 0;
        double donor_angle = 0;
        double system_initial_angle_around_y = 0;
        double acc_disk_solid_color = std::numeric_limits<double>::quiet_NaN();
        double roche_lobe_solid_color = std::numeric_limits<double>::quiet_NaN();
    } config{};

private:
    app() = default;
};

/*
 * following constant is kind of stub, should be removed in future
 *
 * anyway, u're free to increase it
 */
const int MAX_NUMBER_OF_THREADS = 32;


/*
 * binary system main properties and main constants
 */
const double PI = 3.14159265358979323846;

/*
 * distance between stars in R sol
 */
const double L = 0.945;

/*
 * accretor position and properties
 */
const double ACC_X0 = 1;
const double ACC_Y0 = 0;
const double ACC_Z0 = 0;
const double ACC_DISK_R = 0.02;
const double M_ACC = 0.73;

/*
 * donor properties
 */
const double M_DONOR = 0.1;

/*
 * G for Sol mass
 */
const long double G_SOL = 132700000000000000000.;

/*
 * rotational velocity
 */
const double OMEGA = 2 * PI * 10000;
// const long double Msol = 1.988E+30;