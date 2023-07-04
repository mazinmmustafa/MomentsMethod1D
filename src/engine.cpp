//
#include "engine.hpp"

// Engine class

Engine::Engine(){}

Engine::~Engine(){
    Engine::deallocate();
}

void Engine::deallocate(){
    this->Z_mn.deallocate();
    this->V_m.deallocate();
    this->I_n.deallocate();
    this->shape.deallocate();
    this->N_bases = 0;
    this->is_Z_mn_filled = FALSE;
    this->is_V_m_filled = FALSE;
    this->is_I_n_filled = FALSE;
}

// Functions
Basis create_basis(const Vector r_m, const Vector r_n, const Vector r_p){
    Basis basis={r_m, r_n, r_p, Vector(), Vector(), 0.0, 0.0, 0, 0, 0.0, 0.0};
    basis.get_parameters();
    return basis;
}

struct Integral_Args{
    double L=0.0, a=0.0;
    Basis *basis_m=NULL;
    Basis *basis_n=NULL;
    char s_m, s_n;
    Engine *engine=NULL;
    int N_bases=0;
    Vector r;
    Vector direction;
};

// Singular integrals

cmplx I1_integrand_1(cmplx alpha, cmplx alpha_, void *args_in){
    const double k=2.0*pi;
    const cmplx j=cmplx(0.0, 1.0);
    Integral_Args *args=(Integral_Args*)args_in;
    const double L=args->L;
    const double a=args->a;
    const double R=sqrt(pow(L*(real(alpha)-real(alpha_)), 2.0)+a*a);
    const cmplx factor=-j*k*L*L/(4.0*pi);
    return factor*alpha*alpha_*exp(-j*k*R/2.0)*sinc(k*R/2.0);
}

cmplx I1_integrand_2(cmplx alpha, void *args_in){
    Integral_Args *args=(Integral_Args*)args_in;
    const double L=args->L;
    const double a=args->a;
    const cmplx factor=+L/(4.0*pi);
    const double A=sqrt(pow(L*(1.0-real(alpha)), 2.0)+a*a)+L*(1.0-real(alpha));
    const double B=sqrt(pow(L*real(alpha), 2.0)+a*a)-L*real(alpha);
    return factor*pow(alpha, 2.0)*log(A/B);
}

cmplx I1_integrand_3(cmplx alpha, void *args_in){
    Integral_Args *args=(Integral_Args*)args_in;
    const double L=args->L;
    const double a=args->a;
    const cmplx factor=+1.0/(4.0*pi);
    const double A=sqrt(pow(L*(1.0-real(alpha)), 2.0)+a*a);
    const double B=sqrt(pow(L*real(alpha), 2.0)+a*a);
    return factor*real(alpha)*(A-B);
}

cmplx I1_integral(const double L, const double a, QuadL &quadl, int &flag){
    Integral_Args args;
    args.L = L;
    args.a = a;
    cmplx I1=quadl.integral_2D(I1_integrand_1, &args, 0.0, 1.0, 0.0, 1.0, flag);
    cmplx I2=quadl.integral_1D(I1_integrand_2, &args, 0.0, 1.0, flag);
    cmplx I3=quadl.integral_1D(I1_integrand_3, &args, 0.0, 1.0, flag);
    return I1+I2+I3;
}

cmplx I2_integrand_1(cmplx alpha, cmplx alpha_, void *args_in){
    const double k=2.0*pi;
    const cmplx j=cmplx(0.0, 1.0);
    Integral_Args *args=(Integral_Args*)args_in;
    const double L=args->L;
    const double a=args->a;
    const double R=sqrt(pow(L*(real(alpha)-real(alpha_)), 2.0)+a*a);
    const cmplx factor=-j*k/(4.0*pi);
    return factor*exp(-j*k*R/2.0)*sinc(k*R/2.0);
}

cmplx I2_integrand_2(cmplx alpha, void *args_in){
    Integral_Args *args=(Integral_Args*)args_in;
    const double L=args->L;
    const double a=args->a;
    const cmplx factor=+1.0/(4.0*pi*L);
    const double A=sqrt(pow(L*(1.0-real(alpha)), 2.0)+a*a)+L*(1.0-real(alpha));
    const double B=sqrt(pow(L*real(alpha), 2.0)+a*a)-L*real(alpha);
    return factor*log(A/B);
}

cmplx I2_integral(const double L, const double a, QuadL &quadl, int &flag){
    Integral_Args args;
    args.L = L;
    args.a = a;
    cmplx I1=quadl.integral_2D(I2_integrand_1, &args, 0.0, 1.0, 0.0, 1.0, flag);
    cmplx I2=quadl.integral_1D(I2_integrand_2, &args, 0.0, 1.0, flag);
    return I1+I2;
}

cmplx I3_integrand_1(cmplx alpha, cmplx alpha_, void *args_in){
    const double k=2.0*pi;
    const cmplx j=cmplx(0.0, 1.0);
    Integral_Args *args=(Integral_Args*)args_in;
    const double L=args->L;
    const double a=args->a;
    const double R=sqrt(pow(L*(1.0-real(alpha)-real(alpha_)), 2.0)+a*a);
    const cmplx factor=-j*k*L*L/(4.0*pi);
    return factor*alpha*alpha_*exp(-j*k*R/2.0)*sinc(k*R/2.0);
}

cmplx I3_integrand_2(cmplx alpha, void *args_in){
    Integral_Args *args=(Integral_Args*)args_in;
    const double L=args->L;
    const double a=args->a;
    const cmplx factor=+L/(4.0*pi);
    const double A=sqrt(pow(L*(1.0-real(alpha)), 2.0)+a*a)+L*(1.0-real(alpha));
    const double B=sqrt(pow(L*real(alpha), 2.0)+a*a)-L*real(alpha);
    return factor*alpha*(1.0-alpha)*log(A/B);
}

cmplx I3_integrand_3(cmplx alpha, void *args_in){
    Integral_Args *args=(Integral_Args*)args_in;
    const double L=args->L;
    const double a=args->a;
    const cmplx factor=-1.0/(4.0*pi);
    const double A=sqrt(pow(L*(1.0-real(alpha)), 2.0)+a*a);
    const double B=sqrt(pow(L*real(alpha), 2.0)+a*a);
    return factor*real(alpha)*(A-B);
}

cmplx I3_integral(const double L, const double a, QuadL &quadl, int &flag){
    Integral_Args args;
    args.L = L;
    args.a = a;
    cmplx I1=quadl.integral_2D(I3_integrand_1, &args, 0.0, 1.0, 0.0, 1.0, flag);
    cmplx I2=quadl.integral_1D(I3_integrand_2, &args, 0.0, 1.0, flag);
    cmplx I3=quadl.integral_1D(I3_integrand_3, &args, 0.0, 1.0, flag);
    return I1+I2+I3;
}

// Projection
struct Projection_Parameters{
    double l_m=0.0, l_p=0.0;
    double P_m=0.0, P_p=0.0;
    double d=0.0;
};

void projection(Vector v1, Vector v2, Vector p, Projection_Parameters &para){
    check_error(v1.is_equal(v2), "projection is undefined!");
    Vector v21=v2-v1;
    double alpha=real(v21*(p-v1)/pow(v21.mag(), 2.0));
    Vector p0=v1+alpha*v21;
    para.d = (p-p0).mag();
    para.P_m = (p-v1).mag();
    para.P_p = (p-v2).mag();
    if (alpha<0.0){
        para.l_m = +(p0-v1).mag();
        para.l_p = +(p0-v2).mag();
    }else
    if (alpha>1.0){
        para.l_m = -(p0-v1).mag();
        para.l_p = -(p0-v2).mag();
    }else{
        para.l_m = -(p0-v1).mag();
        para.l_p = +(p0-v2).mag();
    }
}

//

double R_mn(const cmplx alpha, const cmplx alpha_, const char s_m, const char s_n, 
    const Basis *basis_m, const Basis *basis_n, double a){
    Vector R;
    if (s_m=='+'&&s_n=='+'){
        R = (basis_m->r_p)-(basis_n->r_p)-alpha*(basis_m->L_p)+alpha_*(basis_n->L_p);
    }else
    if (s_m=='+'&&s_n=='-'){
        R = (basis_m->r_p)-(basis_n->r_m)-alpha*(basis_m->L_p)-alpha_*(basis_n->L_m);
    }else
    if (s_m=='-'&&s_n=='+'){
        R = (basis_m->r_m)-(basis_n->r_p)+alpha*(basis_m->L_m)+alpha_*(basis_n->L_p);
    }else
    if (s_m=='-'&&s_n=='-'){
        R = (basis_m->r_m)-(basis_n->r_m)+alpha*(basis_m->L_m)-alpha_*(basis_n->L_m);
    }else{
        check_error(TRUE, "invalid basis signs option!");
    }
    return sqrt(pow(R.mag(), 2.0)+a*a);
}

cmplx g_mn(const cmplx alpha, const cmplx alpha_, const char s_m, const char s_n, 
    const Basis *basis_m, const Basis *basis_n, double a){
    const double k=2.0*pi;
    const cmplx j=cmplx(0.0, 1.0);
    const double R=R_mn(alpha, alpha_, s_m, s_n, basis_m, basis_n, a);
    return exp(-j*k*R)/(4.0*pi*R);
}

// Case V integrals
cmplx integrand_psi_VI_1(cmplx alpha, cmplx alpha_, void *args_in){
    const double k=2.0*pi;
    const cmplx j=cmplx(0.0, 1.0);
    Integral_Args *args=(Integral_Args*)args_in;
    const Basis *basis_m=args->basis_m;
    const Basis *basis_n=args->basis_n;
    const char s_m=args->s_m;
    const char s_n=args->s_n;
    const double a=args->a;
    cmplx factor=-j*k/(4.0*pi);
    if (s_m=='+'&&s_n=='+'){
        factor*=real(basis_m->L_p*basis_n->L_p);
    }else
    if (s_m=='+'&&s_n=='-'){
        factor*=real(basis_m->L_p*basis_n->L_m);
    }else
    if (s_m=='-'&&s_n=='+'){
        factor*=real(basis_m->L_m*basis_n->L_p);
    }else
    if (s_m=='-'&&s_n=='-'){
        factor*=real(basis_m->L_m*basis_n->L_m);
    }else{
        check_error(TRUE, "invalid basis signs option!");
    }
    double R=R_mn(alpha, alpha_, s_m, s_n, basis_m, basis_n, a);
    return factor*alpha*alpha_*exp(-j*k*R/2.0)*sinc(k*R/2.0);
}

cmplx integrand_psi_VI_2(cmplx alpha, void *args_in){
    Integral_Args *args=(Integral_Args*)args_in;
    const Basis *basis_m=args->basis_m;
    const Basis *basis_n=args->basis_n;
    const char s_m=args->s_m;
    const char s_n=args->s_n;
    const double a=args->a;
    cmplx factor=1.0/(4.0*pi);
    Projection_Parameters para;
    double P_p=0.0, P_m=0.0;
    if (s_m=='+'&&s_n=='+'){
        Vector p=basis_m->r_p-alpha*basis_m->L_p;
        projection(basis_n->r_n, basis_n->r_p, p, para);
        P_p = para.P_p;
        P_m = para.P_m;
        factor*=real(basis_m->L_p*unit(basis_n->L_p))/basis_n->l_p;
        factor*=(-1.0);
    }else
    if (s_m=='+'&&s_n=='-'){
        Vector p=basis_m->r_p-alpha*basis_m->L_p;
        projection(basis_n->r_m, basis_n->r_n, p, para);
        P_p = para.P_p;
        P_m = para.P_m;
        factor*=real(basis_m->L_p*unit(basis_n->L_m))/basis_n->l_m;
        factor*=(+1.0);
    }else
    if (s_m=='-'&&s_n=='+'){
        Vector p=basis_m->r_m+alpha*basis_m->L_m;
        projection(basis_n->r_n, basis_n->r_p, p, para);
        P_p = para.P_p;
        P_m = para.P_m;
        factor*=real(basis_m->L_m*unit(basis_n->L_p))/basis_n->l_p;
        factor*=(-1.0);
    }else
    if (s_m=='-'&&s_n=='-'){
        Vector p=basis_m->r_m+alpha*basis_m->L_m;
        projection(basis_n->r_m, basis_n->r_n, p, para);
        P_p = para.P_p;
        P_m = para.P_m;
        factor*=real(basis_m->L_m*unit(basis_n->L_m))/basis_n->l_m;
        factor*=(+1.0);
    }else{
        check_error(TRUE, "invalid basis signs option!");
    }
    double A=sqrt(pow(P_p, 2.0)+a*a);
    double B=sqrt(pow(P_m, 2.0)+a*a);
    return factor*alpha*(A-B);
}

cmplx integrand_psi_VI_3(cmplx alpha, void *args_in){
    Integral_Args *args=(Integral_Args*)args_in;
    const Basis *basis_m=args->basis_m;
    const Basis *basis_n=args->basis_n;
    const char s_m=args->s_m;
    const char s_n=args->s_n;
    const double a=args->a;
    cmplx factor=1.0/(4.0*pi);
    Projection_Parameters para;
    double P_p=0.0, P_m=0.0;
    double l_p=0.0, l_m=0.0;
    if (s_m=='+'&&s_n=='+'){
        Vector p=basis_m->r_p-alpha*basis_m->L_p;
        projection(basis_n->r_n, basis_n->r_p, p, para);
        P_p = para.P_p;
        P_m = para.P_m;
        l_p = para.l_p;
        l_m = para.l_m;
        factor*=real(basis_m->L_p*unit(basis_n->L_p))/basis_n->l_p;
        factor*=(basis_n->l_p+l_m)/2.0;
    }else
    if (s_m=='+'&&s_n=='-'){
        Vector p=basis_m->r_p-alpha*basis_m->L_p;
        projection(basis_n->r_m, basis_n->r_n, p, para);
        P_p = para.P_p;
        P_m = para.P_m;
        l_p = para.l_p;
        l_m = para.l_m;
        factor*=real(basis_m->L_p*unit(basis_n->L_m))/basis_n->l_m;
        factor*=(-l_m)/2.0;
    }else
    if (s_m=='-'&&s_n=='+'){
        Vector p=basis_m->r_m+alpha*basis_m->L_m;
        projection(basis_n->r_n, basis_n->r_p, p, para);
        P_p = para.P_p;
        P_m = para.P_m;
        l_p = para.l_p;
        l_m = para.l_m;
        factor*=real(basis_m->L_m*unit(basis_n->L_p))/basis_n->l_p;
        factor*=(basis_n->l_p+l_m)/2.0;
    }else
    if (s_m=='-'&&s_n=='-'){
        Vector p=basis_m->r_m+alpha*basis_m->L_m;
        projection(basis_n->r_m, basis_n->r_n, p, para);
        P_p = para.P_p;
        P_m = para.P_m;
        l_p = para.l_p;
        l_m = para.l_m;
        factor*=real(basis_m->L_m*unit(basis_n->L_m))/basis_n->l_m;
        factor*=(-l_m)/2.0;
    }else{
        check_error(TRUE, "invalid basis signs option!");
    }
    double A=(sqrt(pow(P_p, 2.0)+a*a)+l_p)/(sqrt(pow(P_p, 2.0)+a*a)-l_p);
    double B=(sqrt(pow(P_m, 2.0)+a*a)+l_m)/(sqrt(pow(P_m, 2.0)+a*a)-l_m);
    return factor*alpha*(log(A)-log(B));
}

cmplx integrand_phi_VI_1(cmplx alpha, cmplx alpha_, void *args_in){
    const double k=2.0*pi;
    const cmplx j=cmplx(0.0, 1.0);
    Integral_Args *args=(Integral_Args*)args_in;
    const Basis *basis_m=args->basis_m;
    const Basis *basis_n=args->basis_n;
    const char s_m=args->s_m;
    const char s_n=args->s_n;
    const double a=args->a;
    cmplx factor=-j*k/(4.0*pi);
    double R=R_mn(alpha, alpha_, s_m, s_n, basis_m, basis_n, a);
    return factor*exp(-j*k*R/2.0)*sinc(k*R/2.0);
}

cmplx integrand_phi_VI_2(cmplx alpha, void *args_in){
    Integral_Args *args=(Integral_Args*)args_in;
    const Basis *basis_m=args->basis_m;
    const Basis *basis_n=args->basis_n;
    const char s_m=args->s_m;
    const char s_n=args->s_n;
    const double a=args->a;
    cmplx factor=1.0/(4.0*pi);
    Projection_Parameters para;
    double P_p=0.0, P_m=0.0;
    double l_p=0.0, l_m=0.0;
    if (s_m=='+'&&s_n=='+'){
        Vector p=basis_m->r_p-alpha*basis_m->L_p;
        projection(basis_n->r_n, basis_n->r_p, p, para);
        P_p = para.P_p;
        P_m = para.P_m;
        l_p = para.l_p;
        l_m = para.l_m;
        factor*=1.0/basis_n->l_p;
        factor*=1.0/2.0;
    }else
    if (s_m=='+'&&s_n=='-'){
        Vector p=basis_m->r_p-alpha*basis_m->L_p;
        projection(basis_n->r_m, basis_n->r_n, p, para);
        P_p = para.P_p;
        P_m = para.P_m;
        l_p = para.l_p;
        l_m = para.l_m;
        factor*=1.0/basis_n->l_m;
        factor*=1.0/2.0;
    }else
    if (s_m=='-'&&s_n=='+'){
        Vector p=basis_m->r_m+alpha*basis_m->L_m;
        projection(basis_n->r_n, basis_n->r_p, p, para);
        P_p = para.P_p;
        P_m = para.P_m;
        l_p = para.l_p;
        l_m = para.l_m;
        factor*=1.0/basis_n->l_p;
        factor*=1.0/2.0;
    }else
    if (s_m=='-'&&s_n=='-'){
        Vector p=basis_m->r_m+alpha*basis_m->L_m;
        projection(basis_n->r_m, basis_n->r_n, p, para);
        P_p = para.P_p;
        P_m = para.P_m;
        l_p = para.l_p;
        l_m = para.l_m;
        factor*=1.0/basis_n->l_m;
        factor*=1.0/2.0;
    }else{
        check_error(TRUE, "invalid basis signs option!");
    }
    double A=(sqrt(pow(P_p, 2.0)+a*a)+l_p)/(sqrt(pow(P_p, 2.0)+a*a)-l_p);
    double B=(sqrt(pow(P_m, 2.0)+a*a)+l_m)/(sqrt(pow(P_m, 2.0)+a*a)-l_m);
    return factor*(log(A)-log(B));
}

cmplx term_VI(Basis *basis_m, Basis *basis_n, double a, QuadL &quadl, int &flag){
    const double k=2.0*pi;
    const cmplx j=cmplx(0.0, 1.0);
    Integral_Args args;
    args.basis_m = basis_m;
    args.basis_n = basis_n;
    args.a = a;
    cmplx psi_pp, psi_pm, psi_mp, psi_mm;
    cmplx phi_pp, phi_pm, phi_mp, phi_mm;
    args.s_m = '+'; args.s_n='+';
    psi_pp = quadl.integral_2D(integrand_psi_VI_1, &args, 0.0, 1.0, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_psi_VI_2, &args, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_psi_VI_3, &args, 0.0, 1.0, flag);
    phi_pp = quadl.integral_2D(integrand_phi_VI_1, &args, 0.0, 1.0, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_phi_VI_2, &args, 0.0, 1.0, flag);
    args.s_m = '+'; args.s_n='-';
    psi_pm = quadl.integral_2D(integrand_psi_VI_1, &args, 0.0, 1.0, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_psi_VI_2, &args, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_psi_VI_3, &args, 0.0, 1.0, flag);
    phi_pm = quadl.integral_2D(integrand_phi_VI_1, &args, 0.0, 1.0, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_phi_VI_2, &args, 0.0, 1.0, flag);
    args.s_m = '-'; args.s_n='+';
    psi_mp = quadl.integral_2D(integrand_psi_VI_1, &args, 0.0, 1.0, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_psi_VI_2, &args, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_psi_VI_3, &args, 0.0, 1.0, flag);
    phi_mp = quadl.integral_2D(integrand_phi_VI_1, &args, 0.0, 1.0, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_phi_VI_2, &args, 0.0, 1.0, flag);
    args.s_m = '-'; args.s_n='-';
    psi_mm = quadl.integral_2D(integrand_psi_VI_1, &args, 0.0, 1.0, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_psi_VI_2, &args, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_psi_VI_3, &args, 0.0, 1.0, flag);
    phi_mm = quadl.integral_2D(integrand_phi_VI_1, &args, 0.0, 1.0, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_phi_VI_2, &args, 0.0, 1.0, flag);
    cmplx A=psi_pp+psi_pm+psi_mp+psi_mm;
    cmplx B=phi_pp-phi_pm-phi_mp+phi_mm;
    return j*k*eta0*A-j*(eta0/k)*B;
}

// Case VII integrals
cmplx integrand_psi_VII(cmplx alpha, cmplx alpha_, void *args_in){
    Integral_Args *args=(Integral_Args*)args_in;
    const Basis *basis_m=args->basis_m;
    const Basis *basis_n=args->basis_n;
    const char s_m=args->s_m;
    const char s_n=args->s_n;
    const double a=args->a;
    double factor;
    if (s_m=='+'&&s_n=='+'){
        factor = real(basis_m->L_p*basis_n->L_p);
    }else
    if (s_m=='+'&&s_n=='-'){
        factor = real(basis_m->L_p*basis_n->L_m);
    }else
    if (s_m=='-'&&s_n=='+'){
        factor = real(basis_m->L_m*basis_n->L_p);
    }else
    if (s_m=='-'&&s_n=='-'){
        factor = real(basis_m->L_m*basis_n->L_m);
    }else{
        check_error(TRUE, "invalid basis signs option!");
    }
    return factor*alpha*alpha_*g_mn(alpha, alpha_, s_m, s_n, basis_m, basis_n, a);
}

cmplx integrand_phi_VII(cmplx alpha, cmplx alpha_, void *args_in){
    Integral_Args *args=(Integral_Args*)args_in;
    const Basis *basis_m=args->basis_m;
    const Basis *basis_n=args->basis_n;
    const char s_m=args->s_m;
    const char s_n=args->s_n;
    const double a=args->a;
    return g_mn(alpha, alpha_, s_m, s_n, basis_m, basis_n, a);
}

cmplx term_VII(Basis *basis_m, Basis *basis_n, const double a, QuadL &quadl, int &flag){
    const double k=2.0*pi;
    const cmplx j=cmplx(0.0, 1.0);
    Integral_Args args;
    args.basis_m = basis_m;
    args.basis_n = basis_n;
    args.a = a;
    cmplx psi_pp, psi_pm, psi_mp, psi_mm;
    cmplx phi_pp, phi_pm, phi_mp, phi_mm;
    args.s_m = '+'; args.s_n='+';
    psi_pp = quadl.integral_2D(integrand_psi_VII, &args, 0.0, 1.0, 0.0, 1.0, flag);
    phi_pp = quadl.integral_2D(integrand_phi_VII, &args, 0.0, 1.0, 0.0, 1.0, flag);
    args.s_m = '+'; args.s_n='-';
    psi_pm = quadl.integral_2D(integrand_psi_VII, &args, 0.0, 1.0, 0.0, 1.0, flag);
    phi_pm = quadl.integral_2D(integrand_phi_VII, &args, 0.0, 1.0, 0.0, 1.0, flag);
    args.s_m = '-'; args.s_n='+';
    psi_mp = quadl.integral_2D(integrand_psi_VII, &args, 0.0, 1.0, 0.0, 1.0, flag);
    phi_mp = quadl.integral_2D(integrand_phi_VII, &args, 0.0, 1.0, 0.0, 1.0, flag);
    args.s_m = '-'; args.s_n='-';
    psi_mm = quadl.integral_2D(integrand_psi_VII, &args, 0.0, 1.0, 0.0, 1.0, flag);
    phi_mm = quadl.integral_2D(integrand_phi_VII, &args, 0.0, 1.0, 0.0, 1.0, flag);
    cmplx A=psi_pp+psi_pm+psi_mp+psi_mm;
    cmplx B=phi_pp-phi_pm-phi_mp+phi_mm;
    return j*k*eta0*A-j*(eta0/k)*B;
}

// Case I integrals
cmplx term_I(Basis *basis_m, Basis *basis_n, const double a, QuadL &quadl, int &flag){
    const double k=2.0*pi;
    const cmplx j=cmplx(0.0, 1.0);
    Integral_Args args;
    args.basis_m = basis_m;
    args.basis_n = basis_n;
    args.a = a;
    cmplx psi_pp, psi_pm, psi_mp, psi_mm;
    cmplx phi_pp, phi_pm, phi_mp, phi_mm;
    args.s_m = '+'; args.s_n='+';
    psi_pp = I1_integral(basis_m->l_p, a, quadl, flag);
    phi_pp = I2_integral(basis_m->l_p, a, quadl, flag);
    args.s_m = '+'; args.s_n='-';
    psi_pm = quadl.integral_2D(integrand_psi_VI_1, &args, 0.0, 1.0, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_psi_VI_2, &args, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_psi_VI_3, &args, 0.0, 1.0, flag);
    phi_pm = quadl.integral_2D(integrand_phi_VI_1, &args, 0.0, 1.0, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_phi_VI_2, &args, 0.0, 1.0, flag);
    args.s_m = '-'; args.s_n='+';
    psi_mp = quadl.integral_2D(integrand_psi_VI_1, &args, 0.0, 1.0, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_psi_VI_2, &args, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_psi_VI_3, &args, 0.0, 1.0, flag);
    phi_mp = quadl.integral_2D(integrand_phi_VI_1, &args, 0.0, 1.0, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_phi_VI_2, &args, 0.0, 1.0, flag);
    args.s_m = '-'; args.s_n='-';
    psi_mm = I1_integral(basis_m->l_m, a, quadl, flag);
    phi_mm = I2_integral(basis_m->l_m, a, quadl, flag);
    cmplx A=psi_pp+psi_pm+psi_mp+psi_mm;
    cmplx B=phi_pp-phi_pm-phi_mp+phi_mm;
    return j*k*eta0*A-j*(eta0/k)*B;
}

// Case II integrals
cmplx term_II(Basis *basis_m, Basis *basis_n, const double a, QuadL &quadl, int &flag){
    const double k=2.0*pi;
    const cmplx j=cmplx(0.0, 1.0);
    Integral_Args args;
    args.basis_m = basis_m;
    args.basis_n = basis_n;
    args.a = a;
    cmplx psi_pp, psi_pm, psi_mp, psi_mm;
    cmplx phi_pp, phi_pm, phi_mp, phi_mm;
    args.s_m = '+'; args.s_n='+';
    psi_pp = quadl.integral_2D(integrand_psi_VI_1, &args, 0.0, 1.0, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_psi_VI_2, &args, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_psi_VI_3, &args, 0.0, 1.0, flag);
    phi_pp = quadl.integral_2D(integrand_phi_VI_1, &args, 0.0, 1.0, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_phi_VI_2, &args, 0.0, 1.0, flag);
    args.s_m = '+'; args.s_n='-';
    psi_pm = I3_integral(basis_m->l_p, a, quadl, flag);
    phi_pm = I2_integral(basis_m->l_p, a, quadl, flag);
    args.s_m = '-'; args.s_n='+';
    psi_mp = quadl.integral_2D(integrand_psi_VI_1, &args, 0.0, 1.0, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_psi_VI_2, &args, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_psi_VI_3, &args, 0.0, 1.0, flag);
    phi_mp = quadl.integral_2D(integrand_phi_VI_1, &args, 0.0, 1.0, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_phi_VI_2, &args, 0.0, 1.0, flag);
    args.s_m = '-'; args.s_n='-';
    psi_mm = quadl.integral_2D(integrand_psi_VI_1, &args, 0.0, 1.0, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_psi_VI_2, &args, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_psi_VI_3, &args, 0.0, 1.0, flag);
    phi_mm = quadl.integral_2D(integrand_phi_VI_1, &args, 0.0, 1.0, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_phi_VI_2, &args, 0.0, 1.0, flag);
    cmplx A=psi_pp+psi_pm+psi_mp+psi_mm;
    cmplx B=phi_pp-phi_pm-phi_mp+phi_mm;
    return j*k*eta0*A-j*(eta0/k)*B;
}

// Case III integrals
cmplx term_III(Basis *basis_m, Basis *basis_n, const double a, QuadL &quadl, int &flag){
    const double k=2.0*pi;
    const cmplx j=cmplx(0.0, 1.0);
    Integral_Args args;
    args.basis_m = basis_m;
    args.basis_n = basis_n;
    args.a = a;
    cmplx psi_pp, psi_pm, psi_mp, psi_mm;
    cmplx phi_pp, phi_pm, phi_mp, phi_mm;
    args.s_m = '+'; args.s_n='+';
    psi_pp = quadl.integral_2D(integrand_psi_VI_1, &args, 0.0, 1.0, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_psi_VI_2, &args, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_psi_VI_3, &args, 0.0, 1.0, flag);
    phi_pp = quadl.integral_2D(integrand_phi_VI_1, &args, 0.0, 1.0, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_phi_VI_2, &args, 0.0, 1.0, flag);
    args.s_m = '+'; args.s_n='-';
    psi_pm = quadl.integral_2D(integrand_psi_VI_1, &args, 0.0, 1.0, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_psi_VI_2, &args, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_psi_VI_3, &args, 0.0, 1.0, flag);
    phi_pm = quadl.integral_2D(integrand_phi_VI_1, &args, 0.0, 1.0, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_phi_VI_2, &args, 0.0, 1.0, flag);
    args.s_m = '-'; args.s_n='+';
    psi_mp = I3_integral(basis_m->l_m, a, quadl, flag);
    phi_mp = I2_integral(basis_m->l_m, a, quadl, flag);
    args.s_m = '-'; args.s_n='-';
    psi_mm = quadl.integral_2D(integrand_psi_VI_1, &args, 0.0, 1.0, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_psi_VI_2, &args, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_psi_VI_3, &args, 0.0, 1.0, flag);
    phi_mm = quadl.integral_2D(integrand_phi_VI_1, &args, 0.0, 1.0, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_phi_VI_2, &args, 0.0, 1.0, flag);
    cmplx A=psi_pp+psi_pm+psi_mp+psi_mm;
    cmplx B=phi_pp-phi_pm-phi_mp+phi_mm;
    return j*k*eta0*A-j*(eta0/k)*B;
}

// Case IV integrals
cmplx term_IV(Basis *basis_m, Basis *basis_n, const double a, QuadL &quadl, int &flag){
    const double k=2.0*pi;
    const cmplx j=cmplx(0.0, 1.0);
    Integral_Args args;
    args.basis_m = basis_m;
    args.basis_n = basis_n;
    args.a = a;
    cmplx psi_pp, psi_pm, psi_mp, psi_mm;
    cmplx phi_pp, phi_pm, phi_mp, phi_mm;
    args.s_m = '+'; args.s_n='+';
    psi_pp = quadl.integral_2D(integrand_psi_VI_1, &args, 0.0, 1.0, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_psi_VI_2, &args, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_psi_VI_3, &args, 0.0, 1.0, flag);
    phi_pp = quadl.integral_2D(integrand_phi_VI_1, &args, 0.0, 1.0, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_phi_VI_2, &args, 0.0, 1.0, flag);
    args.s_m = '+'; args.s_n='-';
    psi_pm = quadl.integral_2D(integrand_psi_VI_1, &args, 0.0, 1.0, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_psi_VI_2, &args, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_psi_VI_3, &args, 0.0, 1.0, flag);
    phi_pm = quadl.integral_2D(integrand_phi_VI_1, &args, 0.0, 1.0, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_phi_VI_2, &args, 0.0, 1.0, flag);
    args.s_m = '-'; args.s_n='+';
    psi_mp = quadl.integral_2D(integrand_psi_VI_1, &args, 0.0, 1.0, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_psi_VI_2, &args, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_psi_VI_3, &args, 0.0, 1.0, flag);
    phi_mp = quadl.integral_2D(integrand_phi_VI_1, &args, 0.0, 1.0, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_phi_VI_2, &args, 0.0, 1.0, flag);
    args.s_m = '-'; args.s_n='-';
    psi_mm = -I3_integral(basis_m->l_m, a, quadl, flag);
    phi_mm = +I2_integral(basis_m->l_m, a, quadl, flag);
    cmplx A=psi_pp+psi_pm+psi_mp+psi_mm;
    cmplx B=phi_pp-phi_pm-phi_mp+phi_mm;
    return j*k*eta0*A-j*(eta0/k)*B;
}

// Case V integrals
cmplx term_V(Basis *basis_m, Basis *basis_n, const double a, QuadL &quadl, int &flag){
    const double k=2.0*pi;
    const cmplx j=cmplx(0.0, 1.0);
    Integral_Args args;
    args.basis_m = basis_m;
    args.basis_n = basis_n;
    args.a = a;
    cmplx psi_pp, psi_pm, psi_mp, psi_mm;
    cmplx phi_pp, phi_pm, phi_mp, phi_mm;
    args.s_m = '+'; args.s_n='+';
    psi_pp = -I3_integral(basis_m->l_p, a, quadl, flag);
    phi_pp = +I2_integral(basis_m->l_p, a, quadl, flag);
    args.s_m = '+'; args.s_n='-';
    psi_pm = quadl.integral_2D(integrand_psi_VI_1, &args, 0.0, 1.0, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_psi_VI_2, &args, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_psi_VI_3, &args, 0.0, 1.0, flag);
    phi_pm = quadl.integral_2D(integrand_phi_VI_1, &args, 0.0, 1.0, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_phi_VI_2, &args, 0.0, 1.0, flag);
    args.s_m = '-'; args.s_n='+';
    psi_mp = quadl.integral_2D(integrand_psi_VI_1, &args, 0.0, 1.0, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_psi_VI_2, &args, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_psi_VI_3, &args, 0.0, 1.0, flag);
    phi_mp = quadl.integral_2D(integrand_phi_VI_1, &args, 0.0, 1.0, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_phi_VI_2, &args, 0.0, 1.0, flag);
    args.s_m = '-'; args.s_n='-';
    psi_mm = quadl.integral_2D(integrand_psi_VI_1, &args, 0.0, 1.0, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_psi_VI_2, &args, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_psi_VI_3, &args, 0.0, 1.0, flag);
    phi_mm = quadl.integral_2D(integrand_phi_VI_1, &args, 0.0, 1.0, 0.0, 1.0, flag)+
             quadl.integral_1D(integrand_phi_VI_2, &args, 0.0, 1.0, flag);
    cmplx A=psi_pp+psi_pm+psi_mp+psi_mm;
    cmplx B=phi_pp-phi_pm-phi_mp+phi_mm;
    return j*k*eta0*A-j*(eta0/k)*B;
}

int is_near_basis(const Basis *basis_m, const Basis *basis_n, const double near_basis){
    double dist[9];
    dist[0] = (basis_m->r_m-basis_n->r_m).mag();
    dist[1] = (basis_m->r_m-basis_n->r_n).mag();
    dist[2] = (basis_m->r_m-basis_n->r_p).mag();
    dist[3] = (basis_m->r_n-basis_n->r_m).mag();
    dist[4] = (basis_m->r_n-basis_n->r_n).mag();
    dist[5] = (basis_m->r_n-basis_n->r_p).mag();
    dist[6] = (basis_m->r_p-basis_n->r_m).mag();
    dist[7] = (basis_m->r_p-basis_n->r_n).mag();
    dist[8] = (basis_m->r_p-basis_n->r_p).mag();
    double min_dist=dist[0];
    for (int i=1; i<9; i++){
        if (dist[i]<min_dist){
            min_dist = dist[i];
        }
    }
    return min_dist<=near_basis ? TRUE : FALSE;
}

cmplx compute_Z_mn_term(Basis *basis_m, Basis *basis_n, const double a, 
    const double near_basis, QuadL &quadl, int &flag){
    Vector r_m_m=basis_m->r_m;
    Vector r_m_n=basis_m->r_n;
    Vector r_m_p=basis_m->r_p;
    Vector r_n_m=basis_n->r_m;
    Vector r_n_n=basis_n->r_n;
    Vector r_n_p=basis_n->r_p;
    // Case I
    if (r_m_m.is_equal(r_n_m)&&
        r_m_n.is_equal(r_n_n)&&
        r_m_p.is_equal(r_n_p)){
        // std::cout << "I am at case I" << std::endl;
        return term_I(basis_m, basis_n, a, quadl, flag);
    }
    // Case II
    if (r_m_n.is_equal(r_n_m)&&
        r_m_p.is_equal(r_n_n)){
        // std::cout << "I am at case II" << std::endl;
        return term_II(basis_m, basis_n, a, quadl, flag);
    }
    // Case III
    if (r_m_n.is_equal(r_n_p)&&
        r_m_m.is_equal(r_n_n)){
        // std::cout << "I am at case III" << std::endl;
        return term_III(basis_m, basis_n, a, quadl, flag);
    }
    // Case IV
    if (r_m_n.is_equal(r_n_m)&&
        r_m_m.is_equal(r_n_n)){
        // std::cout << "I am at case IV" << std::endl;
        return term_IV(basis_m, basis_n, a, quadl, flag);
    }
    // Case V
    if (r_m_n.is_equal(r_n_p)&&
        r_m_p.is_equal(r_n_n)){
        // std::cout << "I am at case V" << std::endl;
        return term_V(basis_m, basis_n, a, quadl, flag);
    }
    // Case VI
    if (is_near_basis(basis_m, basis_n, near_basis)){
        // std::cout << "I am at case VI" << std::endl;
        return term_VI(basis_m, basis_n, a, quadl, flag);
    }
    // Case VII
    // std::cout << "I am at case VII" << std::endl;
    return term_VII(basis_m, basis_n, a, quadl, flag);
}

// Computer Z_mn

void Engine::fill_Z_mn(){
    if (this->is_Z_mn_filled){
        this->Z_mn.deallocate();
    }
    check_error(!this->shape.get_status(), "shape is not allocated yet!");
    check_error(!this->quadl.get_status(), "quadrature is not set yet!");
    int N_bases;
    double a, lambda;
    this->shape.get_parameters(N_bases, lambda, a);
    this->N_bases = N_bases;
    assert(lambda>0.0);
    double near_basis=this->shape.near_basis;
    this->Z_mn.allocate(N_bases, N_bases);
    Basis *basis_m=NULL, *basis_n=NULL;
    int counter=0;
    int flag = 0;
    for (int m=0; m<N_bases; m++){
        basis_m = this->shape.get_basis(m);
        for (int n=0; n<N_bases; n++){
            basis_n = this->shape.get_basis(n);
            this->Z_mn(m, n) = compute_Z_mn_term(basis_m, basis_n, a, near_basis, 
                this->quadl, flag);
            if (flag){
                printf("\nWarning: element %d, %d did not converge!\n", m, n);
            }
            progress_bar(counter++, N_bases*N_bases, "computing Z_mn matrix...");
        }
    }
    this->is_Z_mn_filled = TRUE;
}

cmplx I_radiation(const cmplx s1, const cmplx s2, const char s_n){
    const cmplx j=cmplx(0.0, 1.0);
    const double k=2.0*pi;
    cmplx A=0.0, B;
    B = pow(k*s1, 2.0);
    if (s_n=='+'){
        A = exp(-j*k*s1)*(1.0+j*k*s1)-1.0;
    }else
    if (s_n=='-'){
        A = exp(+j*k*s1)*(1.0-j*k*s1)-1.0;
    }else{
        check_error(TRUE, "invalid basis signs option!");
    }
    const double tol=1.0E-4;
    return abs(s1)<tol ? exp(+j*k*s2)/2.0 : (A/B)*exp(+j*k*s2);
}

// Compute V_m

void Engine::fill_V_m_delta_gap(){
    check_error(!this->shape.get_status(), "shape is not allocated yet!");
    if (this->is_V_m_filled){
        this->V_m.deallocate();
    }
    if (this->is_I_n_filled){
        this->I_n.deallocate();
    }
    int N_bases;
    double a, lambda;
    this->shape.get_parameters(N_bases, lambda, a);
    assert(lambda*a>0.0);
    this->V_m.allocate(N_bases, 1);
    Basis *basis_m=NULL;
    for (int m=0; m<N_bases; m++){
        basis_m = this->shape.get_basis(m);
        if (basis_m->port_number>0){
            if (basis_m->index==m){
                this->V_m(m, 0) = basis_m->V;
            }
        }else{
            this->V_m(m, 0) = 0.0;
        }
    }
    this->is_V_m_filled = TRUE;
}

void Engine::fill_V_m_plane_wave(const double theta_i, const double phi_i,
        const cmplx E_TM, const cmplx E_TE){
    check_error(!this->shape.get_status(), "shape is not allocated yet!");
    if (this->is_V_m_filled){
        this->V_m.deallocate();
    }
    if (this->is_I_n_filled){
        this->I_n.deallocate();
    }
    int N_bases;
    double a, lambda;
    this->shape.get_parameters(N_bases, lambda, a);
    assert(lambda*a>0.0);
    this->V_m.allocate(N_bases, 1);
    Basis *basis_m=NULL;
    Vector theta_hat(+cos(theta_i)*cos(phi_i),
                     +cos(theta_i)*sin(phi_i),
                     -sin(theta_i));
    Vector phi_hat(-sin(phi_i), +cos(phi_i), +0.0);
    Vector k_hat(sin(theta_i)*cos(phi_i), 
                 sin(theta_i)*sin(phi_i),
                 cos(theta_i));
    Vector L_p, L_m, r_p, r_m;
    cmplx xi_m, xi_p;
    cmplx s1, s2;
    for (int m=0; m<N_bases; m++){
        basis_m = this->shape.get_basis(m);
        L_p = basis_m->L_p;
        L_m = basis_m->L_m;
        r_p = basis_m->r_p;
        r_m = basis_m->r_m;
        s1 = L_p*k_hat;
        s2 = r_p*k_hat;
        xi_p = I_radiation(s1, s2, '+');
        s1 = L_m*k_hat;
        s2 = r_m*k_hat;
        xi_m = I_radiation(s1, s2, '-');
        this->V_m(m, 0) = E_TM*(L_p*theta_hat*xi_p+L_m*theta_hat*xi_m)+
                          E_TE*(L_p*phi_hat*xi_p+L_m*phi_hat*xi_m);
    }
    this->is_V_m_filled = TRUE;
}

// Compute I_n

void Engine::solve_delta_gap(){
    check_error(!this->shape.get_status(), "shape is not allocated yet!");
    int N_bases;
    double a, lambda;
    this->shape.get_parameters(N_bases, lambda, a);
    assert(lambda*a>0.0);
    this->I_n.allocate(N_bases, 1);
    if (!this->is_Z_mn_filled){
        Engine::fill_Z_mn();
    }
    Engine::fill_V_m_delta_gap();
    Basis *basis_m;
    for (int m=0; m<N_bases; m++){
        basis_m = this->shape.get_basis(m);
        if (basis_m->port_number>0){
            if (basis_m->index==m){
                this->Z_mn(m, m)+=basis_m->Z;
            }
        }
    }
    this->I_n.allocate(N_bases, 1);
    this->Z_mn.solve(this->V_m, this->I_n);
    this->is_I_n_filled = TRUE;
    this->fields_quadl.set(this->N_quadl_near_field, this->quadl.get_tol(), this->quadl.get_k_max());
}

void Engine::reset_ports(){
    Basis *basis_m;
    for (int m=0; m<this->N_bases; m++){
        basis_m = this->shape.get_basis(m);
        if (basis_m->port_number>0){
            if (basis_m->index==m){
                this->Z_mn(m, m)-=basis_m->Z;
            }
        }
    }
    this->shape.reset_ports();
}

void Engine::solve_plane_wave(const double theta_i, const double phi_i,
        const cmplx E_TM, const cmplx E_TE){
    check_error(!this->shape.get_status(), "shape is not allocated yet!");
    int N_bases;
    double a, lambda;
    this->shape.get_parameters(N_bases, lambda, a);
    assert(lambda*a>0.0);
    this->I_n.allocate(N_bases, 1);
    if (!this->is_Z_mn_filled){
        Engine::fill_Z_mn();
    }
    Engine::fill_V_m_plane_wave(theta_i, phi_i, E_TM, E_TE);
    this->I_n.allocate(N_bases, 1);
    this->Z_mn.solve(this->V_m, this->I_n);
    this->is_I_n_filled = TRUE;
    this->fields_quadl.set(this->N_quadl_near_field, this->quadl.get_tol(), this->quadl.get_k_max());
}

cmplx Engine::find_Z_in(const int port_number){
    check_error(!this->is_I_n_filled, "there is no solution yet!");
    int index=this->shape.find_port_index(port_number);
    Basis *basis=this->shape.get_basis(index);
    return (basis->V/this->I_n(index, 0))-basis->Z;
}

// Far-Field

Far_Field Engine::get_far_field(const double theta, const double phi){
    check_error(!this->is_I_n_filled, "there is no solution yet!");
    const cmplx j=cmplx(0.0, 1.0);
    const double k=2.0*pi;
    Vector theta_hat(+cos(theta)*cos(phi),
                     +cos(theta)*sin(phi),
                     -sin(theta));
    Vector phi_hat(-sin(phi), +cos(phi), +0.0);
    Vector r_hat(sin(theta)*cos(phi), 
                 sin(theta)*sin(phi),
                 cos(theta));
    int N_bases=this->N_bases;
    cmplx s1, s2;
    Vector sum(0.0, 0.0, 0.0);
    cmplx kappa_p, kappa_m;
    Vector L_p, L_m, r_p, r_m;
    Basis *basis_n;
    for (int n=0; n<N_bases; n++){
        basis_n = this->shape.get_basis(n);
        L_p = basis_n->L_p;
        L_m = basis_n->L_m;
        r_p = basis_n->r_p;
        r_m = basis_n->r_m;
        s1 = L_p*r_hat;
        s2 = r_p*r_hat;
        kappa_p = I_radiation(s1, s2, '+');
        s1 = L_m*r_hat;
        s2 = r_m*r_hat;
        kappa_m = I_radiation(s1, s2, '-');
        sum = sum+this->I_n(n, 0)*(L_p*kappa_p+L_m*kappa_m);
    }
    const double lambda=this->shape.get_lambda();
    Far_Field far_field;
    far_field.E_theta = -j*k*eta0*theta_hat*sum/lambda;
    far_field.E_phi = -j*k*eta0*phi_hat*sum/lambda;
    far_field.H_phi = -j*k*theta_hat*sum/lambda;
    far_field.H_theta = +j*k*phi_hat*sum/lambda;
    return far_field;
}

RCS Engine::get_RCS(const double theta, const double phi){
    check_error(!this->is_I_n_filled, "there is no solution yet!");
    Vector theta_hat(+cos(theta)*cos(phi),
                     +cos(theta)*sin(phi),
                     -sin(theta));
    Vector phi_hat(-sin(phi), +cos(phi), +0.0);
    Vector r_hat(sin(theta)*cos(phi), 
                 sin(theta)*sin(phi),
                 cos(theta));
    int N_bases=this->N_bases;
    cmplx s1, s2;
    cmplx sum_theta=0.0;
    cmplx sum_phi=0.0;
    cmplx kappa_p, kappa_m;
    Vector L_p, L_m, r_p, r_m;
    Basis *basis_n;
    for (int n=0; n<N_bases; n++){
        basis_n = this->shape.get_basis(n);
        L_p = basis_n->L_p;
        L_m = basis_n->L_m;
        r_p = basis_n->r_p;
        r_m = basis_n->r_m;
        s1 = L_p*r_hat;
        s2 = r_p*r_hat;
        kappa_p = I_radiation(s1, s2, '+');
        s1 = L_m*r_hat;
        s2 = r_m*r_hat;
        kappa_m = I_radiation(s1, s2, '-');
        sum_theta+=theta_hat*(kappa_p*L_p+kappa_m*L_m)*this->I_n(n, 0);
        sum_phi+=phi_hat*(kappa_p*L_p+kappa_m*L_m)*this->I_n(n, 0);
    }
    RCS sigma;
    sigma.theta = pi*pow(eta0*abs(sum_theta), 2.0);
    sigma.phi = pi*pow(eta0*abs(sum_phi), 2.0);
    return sigma;
}

// Radiation resistance

cmplx radiation_pattern_U_integrand(const cmplx theta, const cmplx phi, void *args_in){
    Integral_Args *args=(Integral_Args*)args_in;
    Engine *engine=args->engine;
    int N_bases=args->N_bases;
    Vector theta_hat(+cos(theta)*cos(phi),
                     +cos(theta)*sin(phi),
                     -sin(theta));
    Vector phi_hat(-sin(phi), +cos(phi), +0.0);
    Vector r_hat(sin(theta)*cos(phi), 
                 sin(theta)*sin(phi),
                 cos(theta));
    cmplx s1, s2;
    Vector J;
    cmplx kappa_p, kappa_m;
    Vector L_p, L_m, r_p, r_m;
    Basis *basis_n;
    for (int n=0; n<N_bases; n++){
        basis_n =  engine->shape.get_basis(n);
        L_p = basis_n->L_p;
        L_m = basis_n->L_m;
        r_p = basis_n->r_p;
        r_m = basis_n->r_m;
        s1 = L_p*r_hat;
        s2 = r_p*r_hat;
        kappa_p = I_radiation(s1, s2, '+');
        s1 = L_m*r_hat;
        s2 = r_m*r_hat;
        kappa_m = I_radiation(s1, s2, '-');
        J = J+(L_p*kappa_p+L_m*kappa_m)*engine->I_n(n, 0);
    }
    cmplx J_theta=theta_hat*J;
    cmplx J_phi=phi_hat*J;
    return pow(abs(J_theta), 2.0)+pow(abs(J_phi), 2.0);
}

cmplx radiation_pattern_integrand(const cmplx theta, const cmplx phi, void *args_in){
    return radiation_pattern_U_integrand(theta, phi, args_in)*sin(theta);
}

double Engine::radiation_resistance(const int port_number){
    check_error(!this->is_I_n_filled, "there is no solution yet!");
    int index=Engine::shape.find_port_index(port_number);
    cmplx I_in=this->I_n(index, 0);
    Integral_Args args;
    args.engine = this;
    args.N_bases = this->N_bases;
    int flag;
    const double k=2.0*pi;
    double factor=k*k*eta0/(32.0*pi*pi);
    double P_rad=factor*real(this->fields_quadl.integral_2D(radiation_pattern_integrand, &args, 
        0.0, pi, 0.0, 2.0*pi, flag));
    if (flag){printf("Warning: no convergence!\n");}
    return 2.0*P_rad/pow(abs(I_in), 2.0);
}

double Engine::directive_gain(const double theta, const double phi){
    check_error(!this->is_I_n_filled, "there is no solution yet!");
    Integral_Args args;
    args.engine = this;
    args.N_bases = this->N_bases;
    int flag;
    double A=4.0*pi*real(radiation_pattern_U_integrand(theta, phi, &args));
    double B=real(this->fields_quadl.integral_2D(radiation_pattern_integrand, &args, 
        0.0, pi, 0.0, 2.0*pi, flag));
    if (flag){printf("Warning: no convergence!\n");}
    return A/B;
}

// Near field

Vector R_n_vector(const cmplx alpha_, const Vector r, 
        const Basis *basis_n, const char s_n){
    Vector R;
    if (s_n=='+'){
        R = r-basis_n->r_p+alpha_*basis_n->L_p;
    }else
    if (s_n=='-'){
        R = r-basis_n->r_m-alpha_*basis_n->L_m;
    }else{
        check_error(TRUE, "invalid basis signs option!");
    }
    return R;
}

double R_n_mag(const cmplx alpha_, const Vector r, 
        const Basis *basis_n, const char s_n, const double a){
    Vector R=R_n_vector(alpha_, r, basis_n, s_n);
    double R_n_mag_=sqrt(pow(R.mag(), 2.0)+a*a);
    return R_n_mag_;
}

cmplx g_n(const cmplx alpha_, const Vector r, 
        const Basis *basis_n, const char s_n, const double a){
    const cmplx j=cmplx(0.0, 1.0);
    const double k=2.0*pi;
    double R_n=R_n_mag(alpha_, r, basis_n, s_n, a);
    return exp(-j*k*R_n)/(4.0*pi*R_n);
}

cmplx near_field_integral_E(const cmplx alpha_, void *args_in){
    Integral_Args *args=(Integral_Args*)args_in;
    Engine *engine=args->engine;
    int N_bases=args->N_bases;
    Vector r=args->r;
    double a=args->a;
    Vector direction=args->direction;
    const cmplx j=cmplx(0.0, 1.0);
    const double k=2.0*pi;
    Vector I1, I2, I3, I4, I;
    Vector L_p, L_m, r_p, r_m;
    Basis *basis_n;
    double R_n_p, R_n_m;
    Vector R_n_p_hat, R_n_m_hat;
    cmplx A_p, A_m;
    cmplx g_p, g_m;
    for (int n=0; n<N_bases; n++){
        basis_n =  engine->shape.get_basis(n);
        L_p = basis_n->L_p;
        L_m = basis_n->L_m;
        r_p = basis_n->r_p;
        r_m = basis_n->r_m;
        R_n_p = R_n_mag(alpha_, r, basis_n, '+', a);
        R_n_m = R_n_mag(alpha_, r, basis_n, '-', a);
        R_n_p_hat = R_n_vector(alpha_, r, basis_n, '+').unit();
        R_n_m_hat = R_n_vector(alpha_, r, basis_n, '-').unit();
        A_p = (1.0+j*k*R_n_p)/R_n_p;
        A_m = (1.0+j*k*R_n_m)/R_n_m;
        g_p = g_n(alpha_, r, basis_n, '+', a);
        g_m = g_n(alpha_, r, basis_n, '-', a);
        I1 = -j*k*eta0*(L_p*alpha_*g_p);
        I2 = -j*k*eta0*(L_m*alpha_*g_m);
        I3 = -j*(eta0/k)*(R_n_p_hat*A_p*g_p);
        I4 = +j*(eta0/k)*(R_n_m_hat*A_m*g_m);
        I = I+(I1+I2+I3+I4)*engine->I_n(n, 0);
    }
    return I*direction;
}

cmplx near_field_integral_H(const cmplx alpha_, void *args_in){
    Integral_Args *args=(Integral_Args*)args_in;
    Engine *engine=args->engine;
    int N_bases=args->N_bases;
    Vector r=args->r;
    double a=args->a;
    Vector direction=args->direction;
    const cmplx j=cmplx(0.0, 1.0);
    const double k=2.0*pi;
    Vector I1, I2, I;
    Vector L_p, L_m, r_p, r_m;
    Basis *basis_n;
    double R_n_p, R_n_m;
    Vector R_n_p_hat, R_n_m_hat;
    cmplx A_p, A_m;
    cmplx g_p, g_m;
    for (int n=0; n<N_bases; n++){
        basis_n =  engine->shape.get_basis(n);
        L_p = basis_n->L_p;
        L_m = basis_n->L_m;
        r_p = basis_n->r_p;
        r_m = basis_n->r_m;
        R_n_p = R_n_mag(alpha_, r, basis_n, '+', a);
        R_n_m = R_n_mag(alpha_, r, basis_n, '-', a);
        R_n_p_hat = R_n_vector(alpha_, r, basis_n, '+').unit();
        R_n_m_hat = R_n_vector(alpha_, r, basis_n, '-').unit();
        A_p = (1.0+j*k*R_n_p)/R_n_p;
        A_m = (1.0+j*k*R_n_m)/R_n_m;
        g_p = g_n(alpha_, r, basis_n, '+', a);
        g_m = g_n(alpha_, r, basis_n, '-', a);
        I1 = (alpha_*A_p*g_p)*(L_p^R_n_p_hat);
        I2 = (alpha_*A_m*g_m)*(L_m^R_n_m_hat);
        I = I+(I1+I2)*engine->I_n(n, 0);
    }
    return I*direction;
}

Near_Field Engine::get_near_field(const double x, const double y, const double z){
    check_error(!this->is_I_n_filled, "there is no solution yet!");   
    Integral_Args args;
    args.engine = this;
    args.N_bases = this->N_bases;
    int N_bases;
    double lambda;
    double a;
    this->shape.get_parameters(N_bases, lambda, a);
    assert(N_bases*lambda*a>0);
    args.a = a;
    args.r = Vector(x/lambda, y/lambda, z/lambda);
    int flag;
    Vector x_hat(1.0, 0.0, 0.0);
    Vector y_hat(0.0, 1.0, 0.0);
    Vector z_hat(0.0, 0.0, 1.0);
    Near_Field near_field;
    args.direction = x_hat; near_field.E_x = this->fields_quadl.integral_1D(near_field_integral_E, &args, 0.0, 1.0, flag)/lambda;
    // if (flag){printf("Warning: no convergence!\n");}
    args.direction = y_hat; near_field.E_y = this->fields_quadl.integral_1D(near_field_integral_E, &args, 0.0, 1.0, flag)/lambda;
    // if (flag){printf("Warning: no convergence!\n");}
    args.direction = z_hat; near_field.E_z = this->fields_quadl.integral_1D(near_field_integral_E, &args, 0.0, 1.0, flag)/lambda;
    // if (flag){printf("Warning: no convergence!\n");}
    args.direction = x_hat; near_field.H_x = this->fields_quadl.integral_1D(near_field_integral_H, &args, 0.0, 1.0, flag)/lambda;
    // if (flag){printf("Warning: no convergence!\n");}
    args.direction = y_hat; near_field.H_y = this->fields_quadl.integral_1D(near_field_integral_H, &args, 0.0, 1.0, flag)/lambda;
    // if (flag){printf("Warning: no convergence!\n");}
    args.direction = z_hat; near_field.H_z = this->fields_quadl.integral_1D(near_field_integral_H, &args, 0.0, 1.0, flag)/lambda;
    // if (flag){printf("Warning: no convergence!\n");}
    return near_field;
}

// S-parameters

void Engine::get_S_paramters(const cmplx Z_0, Matrix &S){
    int N=this->shape.get_N_ports();
    this->fill_Z_mn();
    for (int m=0; m<N; m++){
        for (int k=0; k<N; k++){
            if (k==m){
                this->shape.add_port(k+1, 1.0, Z_0);
            }else{
                this->shape.add_port(k+1, 0.0, Z_0);
            }
        }
        this->solve_delta_gap();
        for (int n=0; n<N; n++){
            cmplx S_mn;
            if (n==m){
                cmplx Z_in = this->find_Z_in(n+1);
                S(m, n) = (Z_in-Z_0)/(Z_in+Z_0);
            }else{
                int index=this->shape.find_port_index(n+1);
                cmplx I_L=this->I_n(index, 0);
                cmplx V_L=I_L*Z_0; 
                S(m, n) = 2.0*V_L;
            }
        }
        this->reset_ports();
    }    
}
