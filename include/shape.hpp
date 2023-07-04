#ifndef SHAPE_HPP
#define SHAPE_HPP

// Libraries
#include "basic_lib.hpp"
#include "utilities.hpp"
#include "vector.hpp"

using namespace basic_lib;

// Definitions
struct Basis{
    Vector r_m, r_n, r_p;
    Vector L_m, L_p;
    double l_m=0.0, l_p=0.0;
    int index=0;
    int port_number=0;
    cmplx V=0.0, Z=0.0;
    void get_parameters(){
        this->L_m = this->r_n-this->r_m;
        this->L_p = this->r_p-this->r_n;
        this->l_m = this->L_m.mag();
        this->l_p = this->L_p.mag();
    }
};

struct Yag_Antenna{
    int N_reflectors=0;
    int N_directors=0;
    double *reflectors_legnths=NULL;
    double *directors_legnths=NULL;
    double *reflectors_positions=NULL;
    double *directors_positions=NULL;
    double driver_length=0.0;
    double driver_position=0.0;
};

class Shape{
private:
    double lambda=0.0;
    double a=0.0;
    double lc=0.0;
    int N_bases=0;
    int is_allocated=FALSE;
    int is_set=FALSE;
    Basis *bases_list=NULL;
    void allocate();
    int N_ports=0;
public:
    const double near_basis=1.0;
    Shape();
    ~Shape();
    void set_parameters(const double lambda, const double a, const double lc);
    void deallocate();
    void log(const char *filename);
    void get_parameters(int &N_bases, double &lambda, double &a);
    void add_port(const int port_number, const cmplx V, const cmplx Z);
    int find_port_index(const int port_number);
    Basis* get_basis(const int index);
    int get_status(){ return this->is_allocated; }
    int get_N_bases(){ return this->N_bases; }
    double get_lambda(){ return this->lambda; }
    int get_N_ports();
    void reset_ports();
    //
    void create_vertical_dipole(const double L, const int is_port, const int port_number);
    void create_circular_loop(const double r, const int is_port, const int port_number);
    void create_transmission_line(const double L, const double h);
    void create_yagi_antenna(const Yag_Antenna yagi_antnna);
    void create_L_shape(const double L);
    void create_gull_antenna(const double h1, const double h2, 
        const double h3, const double alpha);
    void create_RF_coil(const double R, const double H, const int N_legs);
};

// Functions


#endif