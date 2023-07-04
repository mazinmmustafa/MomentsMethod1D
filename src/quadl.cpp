//
#include "quadl.hpp"

namespace basic_lib{
// Begin library
extern "C" void cgqf_f77_(int *rule, int *order, double *x, double *w);

QuadL::QuadL(){}

QuadL::~QuadL(){
    QuadL::unset();
}

void QuadL::set(const int N_quadl, const double tol, const int k_max){
    check_error(N_quadl<1, "invalid quadrature order!");
    check_error(tol<=0.0, "invalid quadrature tolerance!");
    check_error(k_max<0, "invalid recursion depth!");
    if (this->is_set){
        QuadL::unset();
    }
    this->N_quadl = N_quadl;
    this->tol = tol;
    this->k_max = k_max;
    this->x = (double*)calloc(N_quadl, sizeof(double));
    this->w = (double*)calloc(N_quadl, sizeof(double));
    int rule=1;
    cgqf_f77_(&rule, &(this->N_quadl), this->x, this->w);
    this->is_set = TRUE;
}

void QuadL::unset(){
    if (this->is_set){
        free(this->x);
        free(this->w);
    }
    this->is_set = FALSE;
}

cmplx QuadL::integral_1D_quadl(cmplx func(cmplx, void*), void *args, 
    const double a, const double b){
    double h_p=(b+a)/2.0;
    double h_m=(b-a)/2.0;
    cmplx sum=0.0;
    double w_n;
    double x_n;
    for (int n=0; n<this->N_quadl; n++){
        w_n = this->w[n];
        x_n = this->x[n];
        sum+=h_m*w_n*func(h_m*x_n+h_p, args);
    }
    return sum;
}

cmplx QuadL::integral_1D_recursive(cmplx func(cmplx, void*), void *args, 
    const double a, const double b, int &flag, int k, cmplx I_p){
    double m=(a+b)/2.0;
    cmplx I1=QuadL::integral_1D_quadl(func, args, a, m);
    cmplx I2=QuadL::integral_1D_quadl(func, args, m, b);
    cmplx I_n=I1+I2;
    double error=abs(I_n-I_p);
    double tol=this->tol;
    int k_max=this->k_max;
    if (error>tol*abs(I_n)&&k<this->k_max&&error>0.0){
        k++;
        I1 = QuadL::integral_1D_recursive(func, args, a, m, flag, k, I1);
        I2 = QuadL::integral_1D_recursive(func, args, m, b, flag, k, I2);
        I_n = I1+I2;
    }
    if (k==k_max){
        flag++;
    }
    if (!(flag==0||flag==1)){
        flag = 1;
    }
    return I_n;
}

cmplx QuadL::integral_1D(cmplx func(cmplx, void*), void *args, 
    const double a, const double b, int &flag){
    check_error(!this->is_set, "quadrature is not set yet!");
    flag = 0;
    return QuadL::integral_1D_recursive(func, args, a, b, flag, 0, 0.0);
}

cmplx QuadL::integral_2D_quadl(cmplx func(cmplx, cmplx, void*), void *args, 
    const double a_x, const double b_x, const double a_y, const double b_y){
    double h_x_p=(b_x+a_x)/2.0;
    double h_x_m=(b_x-a_x)/2.0;
    double h_y_p=(b_y+a_y)/2.0;
    double h_y_m=(b_y-a_y)/2.0;
    cmplx sum=0.0;
    double w_m;
    double x_m;
    double w_n;
    double x_n;
    for (int m=0; m<this->N_quadl; m++){
        w_m = this->w[m];
        x_m = this->x[m];
        for (int n=0; n<this->N_quadl; n++){
            w_n = this->w[n];
            x_n = this->x[n];
            sum+=h_x_m*h_y_m*w_m*w_n*func(h_x_m*x_m+h_x_p, h_y_m*x_n+h_y_p, args);
        }
    }
    return sum;
}

cmplx QuadL::integral_2D_recursive(cmplx func(cmplx, cmplx, void*), void *args, 
    const double a_x, const double b_x, const double a_y, const double b_y, int &flag, int k, cmplx I_p){
    double m_x=(a_x+b_x)/2.0;
    double m_y=(a_y+b_y)/2.0;
    cmplx I1=integral_2D_quadl(func, args, a_x, m_x, a_y, m_y);
    cmplx I2=integral_2D_quadl(func, args, a_x, m_x, m_y, b_y);
    cmplx I3=integral_2D_quadl(func, args, m_x, b_x, a_y, m_y);
    cmplx I4=integral_2D_quadl(func, args, m_x, b_x, m_y, b_y);
    cmplx I_n=I1+I2+I3+I4;
    double error=abs(I_n-I_p);
    double tol=this->tol;
    int k_max=this->k_max;
    if (error>tol*abs(I_n)&&k<k_max&&error>0.0){
        k++;
        I1 = integral_2D_recursive(func, args, a_x, m_x, a_y, m_y, flag, k, I1);
        I2 = integral_2D_recursive(func, args, a_x, m_x, m_y, b_y, flag, k, I2);
        I3 = integral_2D_recursive(func, args, m_x, b_x, a_y, m_y, flag, k, I3);
        I4 = integral_2D_recursive(func, args, m_x, b_x, m_y, b_y, flag, k, I4);
        I_n = I1+I2+I3+I4; 
    }
    if (k==k_max){
        flag++;
    }
    if (!(flag==0||flag==1)){
        flag = 1;
    }
    return I_n;
}

cmplx QuadL::integral_2D(cmplx func(cmplx, cmplx, void*), void *args, 
    const double a_x, const double b_x, const double a_y, const double b_y, int &flag){
    check_error(!this->is_set, "quadrature is not set yet!");
    flag = 0;
    return QuadL::integral_2D_recursive(func, args, a_x, b_x, a_y, b_y, flag, 0, 0.0);
}

// End of library
}