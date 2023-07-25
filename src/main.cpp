// 
#include "testbench.hpp"
#include "programs.hpp"

int main(){

    // Testbench
    test_utilities();
    test_range();
    test_bessel();
    test_matrix();
    test_vector();
    test_quadl();
    test_shape();
    test_integrals();
    test_engine();

    // Programs
    // vertical_dipole_input_admittance();
    // vertical_dipole_radiation_resistance();
    // vertical_dipole_directive_gain();
    // vertical_dipole_near_field_2D();
    // vertical_dipole_near_field_vs_far_field();
    // vertical_dipole_RCS_1();
    // vertical_dipole_RCS_2();
    // vertical_dipole_RCS_3();
    // monopole_s_parameters();
    // L_shape_RCS_1();
    // L_shape_RCS_2();
    // circular_loop_input_impedance();
    // transmission_line_input_impedance();
    // transmission_line_s_parameters();
    // transmission_line_near_field();
    // transmission_line_far_field();
    // yagi_uda_antenna_far_field(); 
    // gull_antenna_far_field();
    // RF_coil_s_parameters();

    return 0;
}
