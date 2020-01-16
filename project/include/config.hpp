#pragma once

const int AMOUNT_OF_THREADS = 16;
const double PI = 3.14159265358979323846;
const double marked_accretor_value = std::numeric_limits<double>::quiet_NaN();

/*
 * TODO:
 * 1. move consts to singleton
 */
const double L = 0.945;

// accretor properties
const double acc_x0 = 1;
const double acc_y0 = 0;
const double acc_z0 = 0;

const double acc_r = 0.02;
const double Ma = 0.73;

// donor properties
const double Md = 0.1;

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