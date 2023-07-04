#ifndef ENGINE_HPP
#define ENGINE_HPP

// Libraries
#include "basic_lib.hpp"
#include "utilities.hpp"
#include "quadl.hpp"
#include "matrix.hpp"
#include "vector.hpp"
#include "shape.hpp"

using namespace basic_lib;

// Definitions
struct Far_Field{
    cmplx E_theta=0.0, E_phi=0.0;
    cmplx H_theta=0.0, H_phi=0.0;
};

struct Near_Field{
    cmplx E_x=0.0, E_y=0.0, E_z=0.0;
    cmplx H_x=0.0, H_y=0.0, H_z=0.0;
};

struct RCS{
    double theta=0.0, phi=0.0;
};

class Engine{
private:
    int N_bases=0;
    int is_Z_mn_filled=FALSE;
    int is_V_m_filled=FALSE;
    int is_I_n_filled=FALSE;
    QuadL fields_quadl;
    const int N_quadl_near_field=32;
public:
    Engine();
    ~Engine();
    Matrix Z_mn, V_m, I_n;
    QuadL quadl;
    Shape shape;
    void deallocate();
    void fill_Z_mn();
    void fill_V_m_delta_gap();
    void reset_ports();
    void fill_V_m_plane_wave(const double theta_i, const double phi_i,
        const cmplx E_TM, const cmplx E_TE);
    void solve_delta_gap();
    void solve_plane_wave(const double theta_i, const double phi_i,
        const cmplx E_TM, const cmplx E_TE);
    cmplx find_Z_in(const int port_number);
    Far_Field get_far_field(const double theta, const double phi);
    RCS get_RCS(const double theta, const double phi);
    double radiation_resistance(const int port_number);
    double directive_gain(const double theta, const double phi);
    Near_Field get_near_field(const double x, const double y, const double z);
    void get_S_paramters(const cmplx Z_0, Matrix &S);
};

// Functions
Basis create_basis(const Vector r_m, const Vector r_n, const Vector r_p);
cmplx I1_integral(const double L, const double a, QuadL &quadl, int &flag);
cmplx I2_integral(const double L, const double a, QuadL &quadl, int &flag);
cmplx I3_integral(const double L, const double a, QuadL &quadl, int &flag);
cmplx term_I(Basis *basis_m, Basis *basis_n, const double a, QuadL &quadl, int &flag);
cmplx term_II(Basis *basis_m, Basis *basis_n, const double a, QuadL &quadl, int &flag);
cmplx term_III(Basis *basis_m, Basis *basis_n, const double a, QuadL &quadl, int &flag);
cmplx term_IV(Basis *basis_m, Basis *basis_n, const double a, QuadL &quadl, int &flag);
cmplx term_V(Basis *basis_m, Basis *basis_n, const double a, QuadL &quadl, int &flag);
cmplx term_VI(Basis *basis_m, Basis *basis_n, const double a, QuadL &quadl, int &flag);
cmplx term_VII(Basis *basis_m, Basis *basis_n, const double a, QuadL &quadl, int &flag);
cmplx compute_Z_mn_term(Basis *basis_m, Basis *basis_n, const double a, 
    const double near_basis, QuadL &quadl, int &flag);

#endif