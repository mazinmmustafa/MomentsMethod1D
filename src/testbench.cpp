//
#include "testbench.hpp"

using namespace basic_lib;

void test_utilities(){

    printf("Testing utilities:\n");

    check_error(1==0, "This is not an error!");

    File file;

    file.open("data/misc/test_file.txt", "w");
    file.write("%21.14E %21.14E %21.14E\n", pi, c0, mu0);
    file.write("%21.14E %21.14E %21.14E\n", eps0, pi, eta0);
    file.close();

    file.open("data/misc/test_file.txt", "a");
    file.write("The next line is appended:\n");
    file.write("%21.14E %21.14E %21.14E\n", pi, c0, mu0);
    file.close();

    file.open("data/misc/test_file.txt", "r");
    double x, y;
    file.read("%lf %lf", &x, &y);
    printf("%21.14E %21.14E\n", x, y);
    file.read("%lf", &x);
    printf("%21.14E\n", x);
    file.read("%lf", &x);
    printf("%21.14E\n", x);
    file.close();

    Timer T;
    T.unset();
    T.set();
    
    int N=600;
    Random random_gen;
    double r_min=-4.0, r_max=+3.0;

    File new_file;
    new_file.open("data/misc/randon.dat", "w");

    T.set();
    for (int i=0; i<N; i++){
        new_file.write("%3d %21.14E\n", i, random_gen.rand_double(r_min, r_max));
        usleep(1000);
        progress_bar(i, N, " dummy loop");
    }
    T.unset();

    new_file.close();

}

void test_range(){
    
    printf("Testing range:\n");

    Range x_range;
    File file;

    x_range.linspace(0.0, 10.0, 60);
    file.open("data/test_basic_lib/test_range/data1.dat", "w");
    for (int i=0; i<x_range.size(); i++){
        double x, y;
        x = x_range(i);
        y = sin(sin(x)+x);
        file.write("%21.14E %21.14E\n", x, y);
    }
    file.close();

    x_range.logspace(1.0E-3, 1.0E+1, 50);
    file.open("data/test_basic_lib/test_range/data2.dat", "w");
    for (int i=0; i<x_range.size(); i++){
        double x, y;
        x = x_range(i);
        y = sin(x)*exp(-x);
        file.write("%21.14E %21.14E\n", x, y);
    }
    file.close();

}

void test_bessel(){

    printf("Testing bessel:\n");

    int Ns=100;
    double R=4.0;

    Range x_range, y_range;
    x_range.linspace(-R, +R, Ns);
    y_range.linspace(-R, +R, Ns);

    File file;
    file.open("data/test_basic_lib/test_bessel/data1.dat", "w");
    double x, y;
    cmplx z;
    for (int i=0; i<x_range.size(); i++){
        x = x_range(i);
        for (int j=0; j<y_range.size(); j++){
            y = y_range(j);
            z = cmplx(x, y);
            file.write("%21.14E %21.14E %21.14E\n", x, y, log10(abs(besselh(2, 1, z))));
        }
        file.write("\n");
    }
    file.close();

}

void test_matrix(){

    printf("Testing matrix:\n");

    Matrix A;
    A.allocate(4, 3);
    A.eye();

    A.save("data/misc/A.dat");

    Matrix B;
    B.allocate(4, 3);
    B.load("data/misc/A.dat");

    disp(A(0, 0));
    B(1, 0) = -1.0;
    disp(B(1, 0));
    A.deallocate();

    A.allocate(3, 3);
    A(0, 0) = +0.0; A(0, 1) = +1.0; A(0, 2) = +2.0; 
    A(1, 0) = -4.0; A(1, 1) = +0.0; A(1, 2) = +0.0; 
    A(2, 0) = +6.0; A(2, 1) = +0.0; A(2, 2) = +1.0; 

    A.lup();
    A.lup();
    B.allocate(3, 3);
    A.copy_lu(B);
    B.save("data/misc/LU.dat");

    Matrix b, x;
    b.allocate(3, 1);
    x.allocate(3, 1);
    b(0, 0) = +1.0; 
    b(1, 0) = +2.0;
    b(2, 0) = +3.0;

    A.solve(b, x);
    x.save("data/misc/x.dat");
    disp(A.det());
    A.inv();
    A.save("data/misc/inv.dat");
    disp(A.det());

}

void test_vector(){

    printf("Testing vector:\n");

    Vector A(1.0, 2.0, 3.0);
    Vector B=A;
    disp(A^B);
    disp(A+B);
    disp(A-B);
    disp(A+4.0*B/3.0-0.5*A.unit());
    disp(A*B);
    disp(mag(A)*mag(A));
    disp(pow(A*unit(A), 2.0));
    
}

struct Func_Args{
    double a;
};

cmplx func_1D(const cmplx x, void *args){
    Func_Args *Args=(Func_Args*)args;
    double a=Args->a;
    return exp(-a*a*x*x)+sin(x);
}

cmplx func_2D(const cmplx x, const cmplx y, void *args){
    Func_Args *Args=(Func_Args*)args;
    double a=Args->a;
    return cos(x+1.0)*log(2.0-y)+exp(-a*a*(pow(x-0.2, 2.0)+pow(y-0.2, 2.0)));
}

void test_quadl(){ 

    printf("Testing quadl:\n");

    int N_quadl=4;
    double tol=1.0E-12;
    int k_max=15;
    int flag;
    QuadL quadl;
    quadl.set(N_quadl, tol, k_max);
    quadl.set(N_quadl, tol, k_max);

    Func_Args args;
    Timer T;
    
    args.a = 10.0;
    T.set();
    disp(quadl.integral_1D(func_1D, &args, -2.0, +4.0, flag));
    T.unset();
    printf("flag = %d\n", flag);

    args.a = sqrt(80.0);
    T.set();
    disp(quadl.integral_2D(func_2D, &args, -2.0, 1.0, -1.0, 1.0, flag));
    T.unset();
    printf("flag = %d\n", flag);

    quadl.unset();
    quadl.unset();

}

void test_shape(){

    double lambda=1.0;

    double L=0.5*lambda;
    double a=1.0E-3*lambda;
    double lc=1.0/21.0;

    Shape shape;
    shape.set_parameters(lambda, a, lc);
    shape.create_vertical_dipole(L, TRUE, 1);

    std::cout << shape.find_port_index(1) << std::endl;
    shape.add_port(1, 1.0, 50.0);

    shape.log("data/misc/shape_log.txt");  
    
}

void test_integrals(){

    int N_quadl=4;
    double tol=1.0E-6;
    int k_max=15;
    QuadL quadl;
    quadl.set(N_quadl, tol, k_max);
    int flag;

    double lambda=1.0;
    double L=0.5*lambda;
    double a=1.0E-5*lambda;
    double lc=1.0/21.0;

    Engine engine;
    engine.shape.set_parameters(lambda, a, lc);
    engine.shape.create_vertical_dipole(L, TRUE, 1);
    engine.shape.log("data/misc/shape_log.txt");

    L = 0.1;
    a = 1.0E-4;
    std::cout << "Singular integral 1:" << std::endl;
    disp(I1_integral(L, a, quadl, flag)); disp(flag);

    std::cout << "Singular integral 2:" << std::endl;
    disp(I2_integral(L, a, quadl, flag)); disp(flag);

    std::cout << "Singular integral 3:" << std::endl;
    disp(I3_integral(L, a, quadl, flag)); disp(flag);

    L = 1.0/15.0;
    a = 1.0E-5;

    Vector r_m_m, r_m_n, r_m_p;
    Vector r_n_m, r_n_n, r_n_p;

    Basis basis_m;
    Basis basis_n;

    cmplx ans;

    Timer T;
    T.set();

    // Case I
    std::cout << "Case I:" << std::endl;

    r_m_m.create_vector(+0.0*L, +0.0*L, +0.0*L);
    r_m_n.create_vector(+0.2*L, +0.0*L, +0.0*L);
    r_m_p.create_vector(+0.6*L, +0.0*L, +0.0*L);

    r_n_m.create_vector(+0.0*L, +0.0*L, +0.0*L);
    r_n_n.create_vector(+0.2*L, +0.0*L, +0.0*L);
    r_n_p.create_vector(+0.6*L, +0.0*L, +0.0*L);

    basis_m = create_basis(r_m_m, r_m_n, r_m_p);
    basis_n = create_basis(r_n_m, r_n_n, r_n_p);
    
    ans = term_I(&basis_m, &basis_n, a, quadl, flag);
    disp(ans);
    disp(flag);
    // disp(term_VII(&basis_m, &basis_n, a, quadl, flag));
    // disp(flag);
    disp(compute_Z_mn_term(&basis_m, &basis_n, a, 1.0, quadl, flag));

    // Case II
    std::cout << "Case II:" << std::endl;

    r_m_m.create_vector(+0.0*L, +0.0*L, +0.0*L);
    r_m_n.create_vector(+0.2*L, +0.0*L, +0.0*L);
    r_m_p.create_vector(+0.6*L, +0.0*L, +0.0*L);

    r_n_m.create_vector(+0.2*L, +0.0*L, +0.0*L);
    r_n_n.create_vector(+0.6*L, +0.0*L, +0.0*L);
    r_n_p.create_vector(+0.12*L, +0.0*L, +0.0*L);

    basis_m = create_basis(r_m_m, r_m_n, r_m_p);
    basis_n = create_basis(r_n_m, r_n_n, r_n_p);
    
    ans = term_II(&basis_m, &basis_n, a, quadl, flag);
    disp(ans);
    disp(flag);
    // disp(term_VII(&basis_m, &basis_n, a, quadl, flag));
    // disp(flag);
    disp(compute_Z_mn_term(&basis_m, &basis_n, a, 1.0, quadl, flag));

    // Case III
    std::cout << "Case III:" << std::endl;

    r_m_m.create_vector(+0.2*L, +0.0*L, +0.0*L);
    r_m_n.create_vector(+0.6*L, +0.0*L, +0.0*L);
    r_m_p.create_vector(+0.12*L, +0.0*L, +0.0*L);

    r_n_m.create_vector(+0.0*L, +0.0*L, +0.0*L);
    r_n_n.create_vector(+0.2*L, +0.0*L, +0.0*L);
    r_n_p.create_vector(+0.6*L, +0.0*L, +0.0*L);

    basis_m = create_basis(r_m_m, r_m_n, r_m_p);
    basis_n = create_basis(r_n_m, r_n_n, r_n_p);
    
    ans = term_III(&basis_m, &basis_n, a, quadl, flag);
    disp(ans);
    disp(flag);
    // disp(term_VII(&basis_m, &basis_n, a, quadl, flag));
    // disp(flag);
    disp(compute_Z_mn_term(&basis_m, &basis_n, a, 1.0, quadl, flag));

    // Case IV
    std::cout << "Case IV:" << std::endl;

    r_m_m.create_vector(+0.0*L, +0.0*L, +0.0*L);
    r_m_n.create_vector(+0.3*L, +0.0*L, +0.0*L);
    r_m_p.create_vector(+0.7*L, +0.0*L, +0.0*L);

    r_n_m.create_vector(+0.3*L, +0.0*L, +0.0*L);
    r_n_n.create_vector(+0.0*L, +0.0*L, +0.0*L);
    r_n_p.create_vector(-0.6*L, +0.0*L, +0.0*L);

    basis_m = create_basis(r_m_m, r_m_n, r_m_p);
    basis_n = create_basis(r_n_m, r_n_n, r_n_p);
    
    ans = term_IV(&basis_m, &basis_n, a, quadl, flag);
    disp(ans);
    disp(flag);
    // disp(term_VII(&basis_m, &basis_n, a, quadl, flag));
    // disp(flag);
    disp(compute_Z_mn_term(&basis_m, &basis_n, a, 1.0, quadl, flag));

    // Case V
    std::cout << "Case V:" << std::endl;

    r_m_m.create_vector(+0.7*L, +0.0*L, +0.0*L);
    r_m_n.create_vector(+0.3*L, +0.0*L, +0.0*L);
    r_m_p.create_vector(+0.0*L, +0.0*L, +0.0*L);

    r_n_m.create_vector(-0.6*L, +0.0*L, +0.0*L);
    r_n_n.create_vector(+0.0*L, +0.0*L, +0.0*L);
    r_n_p.create_vector(+0.3*L, +0.0*L, +0.0*L);

    basis_m = create_basis(r_m_m, r_m_n, r_m_p);
    basis_n = create_basis(r_n_m, r_n_n, r_n_p);
    
    ans = term_V(&basis_m, &basis_n, a, quadl, flag);
    disp(ans);
    disp(flag);
    // disp(term_VII(&basis_m, &basis_n, a, quadl, flag));
    // disp(flag);
    disp(compute_Z_mn_term(&basis_m, &basis_n, a, 1.0, quadl, flag));

    // Case VI
    std::cout << "Case VI:" << std::endl;

    r_m_m.create_vector(-0.8*L, +0.0*L, +0.0*L);
    r_m_n.create_vector(-0.3*L, +0.0*L, +0.0*L);
    r_m_p.create_vector(+0.0*L, +0.0*L, +0.0*L);

    r_n_m.create_vector(+0.0*L, +0.0*L, +0.0*L);
    r_n_n.create_vector(+0.6*L, +0.0*L, +0.0*L);
    r_n_p.create_vector(+0.9*L, +0.0*L, +0.0*L);

    basis_m = create_basis(r_m_m, r_m_n, r_m_p);
    basis_n = create_basis(r_n_m, r_n_n, r_n_p);
    
    ans = term_VI(&basis_m, &basis_n, a, quadl, flag);
    disp(ans);
    disp(flag);
    disp(compute_Z_mn_term(&basis_m, &basis_n, a, 1.0, quadl, flag));

    // Case VII
    std::cout << "Case VII:" << std::endl;

    r_m_m.create_vector(-0.8*L, +0.0*L, +0.0*L);
    r_m_n.create_vector(-0.3*L, +0.0*L, +0.0*L);
    r_m_p.create_vector(+0.0*L, +0.0*L, +0.0*L);

    r_n_m.create_vector(+0.0*L, +0.0*L, +2.0*lambda);
    r_n_n.create_vector(+0.6*L, +0.0*L, +2.0*lambda);
    r_n_p.create_vector(+0.9*L, +0.0*L, +2.0*lambda);

    basis_m = create_basis(r_m_m, r_m_n, r_m_p);
    basis_n = create_basis(r_n_m, r_n_n, r_n_p);
    
    ans = term_VII(&basis_m, &basis_n, a, quadl, flag);
    disp(ans);
    disp(flag);
    disp(compute_Z_mn_term(&basis_m, &basis_n, a, 1.0, quadl, flag));
    
    T.unset();

}

void test_engine(){

    // Quadrature
    int N_quadl=4;
    double tol=1.0E-4;
    int k_max=15;

    // Definitions
    double lambda=1.0;
    double L=0.47*lambda;
    double a=5.0E-3*lambda;
    cmplx V_in=1.0;
    cmplx Z_0=50.0;
    double lc=1.0/21.0;

    // Engine
    Timer T;
    Engine engine;
    engine.quadl.set(N_quadl, tol, k_max);
    engine.shape.set_parameters(lambda, a, lc);
    engine.shape.create_vertical_dipole(L, TRUE, 1);
    engine.shape.add_port(1, V_in, Z_0);
    engine.shape.log("data/misc/shape_log.txt");
    T.set();
    engine.solve_delta_gap();
    T.unset();
    engine.Z_mn.save("data/misc/Z_mn.dat");
    engine.V_m.save("data/misc/V_m.dat");
    engine.I_n.save("data/misc/I_n.dat");

    File file;
    file.open("data/misc/I_dipole.dat", "w");
    int N_bases=engine.shape.get_N_bases();
    Basis *basis_n;
    for (int n=0; n<N_bases; n++){
        basis_n=engine.shape.get_basis(n);
        file.write("%21.14E %21.14E %21.14E\n", real(basis_n->r_n.z), 
            abs(engine.I_n(n, 0))/1.0E-3, rad2deg(arg(engine.I_n(n, 0))));
    }
    file.close();

    disp(engine.find_Z_in(1));

}