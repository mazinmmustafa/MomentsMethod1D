#ifndef QUADL_HPP
#define QUADL_HPP

// Libraries
#include "basic_lib.hpp"
#include "utilities.hpp"

namespace basic_lib{
// Definitions
class QuadL{
private:
    int N_quadl=0;
    int k_max=0;
    double tol=0.0;
    double *x=NULL;
    double *w=NULL;
    int is_set=FALSE;
    //
    cmplx integral_1D_quadl(cmplx func(cmplx, void*), 
        void *args, const double a, const double b);
    cmplx integral_1D_recursive(cmplx func(cmplx, void*), 
        void *args, const double a, double b, int &flag, int k, cmplx I_p);
    cmplx integral_2D_quadl(cmplx func(cmplx, cmplx, void*), 
        void *args, const double a_x, const double b_x, const double a_y, const double b_y);
    cmplx integral_2D_recursive(cmplx func(cmplx, cmplx, void*), 
        void *args, const double a_x, const double b_x, const double a_y, const double b_y, 
        int &flag, int k, cmplx I_p);
public:
    QuadL();
    ~QuadL();
    int get_status(){return this->is_set;}
    double get_tol(){return this->tol;}
    int get_k_max(){return this->k_max;}
    void set(const int N_quadl, const double tol, const int k_max);
    void unset();
    cmplx integral_1D(cmplx func(cmplx, void*), void *args, const double a, const double b, int &flag);
    cmplx integral_2D(cmplx func(cmplx, cmplx, void*), void *args, 
        const double a_x, const double b_x, const double a_y, const double b_y, int &flag);
};

// Functions

}

#endif