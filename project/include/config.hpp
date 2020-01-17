#pragma once

const int AMOUNT_OF_THREADS = 16;
const double PI = 3.14159265358979323846;
const double ACC_DISK_COLOR = std::numeric_limits<double>::quiet_NaN();
const double ROCHE_LOBE_COLOR = std::numeric_limits<double>::quiet_NaN();

/*
 * TODO:
 * 1. move consts to singleton
 */
const double L = 0.945;

// accretor properties
const double ACC_X0 = 1;
const double ACC_Y0 = 0;
const double ACC_Z0 = 0;

const double ACC_DISK_R = 0.02;
const double M_ACC = 0.73;

// donor properties
const double M_DONOR = 0.1;

/*
 * for Sol mass
 */
const long double G_SOL = 132700000000000000000.;
const double OMEGA = 2 * PI * 10000;
// const long double Msol = 1.988E+30;

/*
 * used to select donor position
 */
const double angle_around_x_axis = 0;
const double angle_around_y_axis = 0;