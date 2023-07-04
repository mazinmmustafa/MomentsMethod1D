#ifndef PROGRAMS_HPP
#define PROGRAMS_HPP

// Libraries
#include "basic_lib.hpp"
#include "utilities.hpp"
#include "bessel.hpp"
#include "matrix.hpp"
#include "vector.hpp"
#include "quadl.hpp"
#include "shape.hpp"
#include "engine.hpp"

// Definitions

// Functions
void vertical_dipole_input_admittance();
void circular_loop_input_impedance();
void transmission_line_input_impedance();
void transmission_line_s_parameters();
void yagi_uda_antenna_far_field();
void vertical_dipole_RCS_1();
void vertical_dipole_RCS_2();
void vertical_dipole_RCS_3();
void monopole_s_parameters();
void L_shape_RCS_1();
void L_shape_RCS_2();
void gull_antenna_far_field();
void vertical_dipole_radiation_resistance();
void vertical_dipole_directive_gain();
void transmission_line_near_field();
void vertical_dipole_near_field_2D();
void vertical_dipole_near_field_vs_far_field();
void transmission_line_far_field();
void RF_coil_s_parameters();

#endif