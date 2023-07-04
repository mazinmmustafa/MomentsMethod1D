//
#include "programs.hpp"

using namespace basic_lib;

void vertical_dipole_input_admittance(){

    Timer T;
    T.set();

    // Units
    const double mS=1.0E-3;

    // Definitions
    const int Ns=801;
    const int N_quadl=2;
    const double tol=1.0E-3;
    const int k_max=15;
    const double lambda=1.0;
    const double a=1.0E-3*lambda;
    const double lc=1.0/21.0;
    const cmplx V_in=1.0;
    const cmplx Z_0=50.0;
    const double L_min=0.01*lambda;
    const double L_max=4.00*lambda;

    // Program
    Range L;    
    L.linspace(L_min, L_max, Ns);
    File file;
    Engine engine;
    engine.quadl.set(N_quadl, tol, k_max);
    file.open("data/vertical_dipole/input_admittance/data1.dat", "w");
    cmplx Z_in, Y_in;
    for (int i=0; i<Ns; i++){
        std::cout << "Step " << i+1 << " out of " << Ns << std::endl;
        engine.shape.set_parameters(lambda, a, lc);
        engine.shape.create_vertical_dipole(L(i), TRUE, 1);
        engine.shape.add_port(1, V_in, Z_0);
        engine.solve_delta_gap();
        Z_in = engine.find_Z_in(1);
        Y_in = 1.0/Z_in;
        file.write("%21.14E %21.14E %21.14E\n", L(i)/lambda, real(Y_in)/mS, imag(Y_in)/mS);
        engine.deallocate();
    }
    file.close();

    T.unset();

}

void circular_loop_input_impedance(){

    Timer T;
    T.set();

    // Units
    const double k_ohm=1.0E+3;

    // Definitions
    const int Ns=801;
    const int N_quadl=2;
    const double tol=1.0E-3;
    const int k_max=15;
    const double lambda=1.0;
    const double a=1.0E-4*lambda;
    const double lc=1.0/21.0;
    const cmplx V_in=1.0;
    const cmplx Z_0=50.0;
    const double S_min=0.01*lambda;
    const double S_max=4.00*lambda;

    // Program
    Range S;    
    S.linspace(S_min, S_max, Ns);
    File file;
    Engine engine;
    engine.quadl.set(N_quadl, tol, k_max);
    file.open("data/circular_loop/input_impedance/data1.dat", "w");
    cmplx Z_in;
    for (int i=0; i<Ns; i++){
        std::cout << "Step " << i+1 << " out of " << Ns << std::endl;
        engine.shape.set_parameters(lambda, a, lc);
        engine.shape.create_circular_loop(S(i)/(2.0*pi), TRUE, 1);
        engine.shape.add_port(1, V_in, Z_0);
        engine.solve_delta_gap();
        Z_in = engine.find_Z_in(1);
        file.write("%21.14E %21.14E %21.14E\n", S(i)/lambda, real(Z_in)/k_ohm, imag(Z_in)/k_ohm);
        engine.deallocate();
    }
    file.close();

    T.unset();

}

void transmission_line_input_impedance(){

    Timer T;
    T.set();

    // Units
    const double kHz=1.0E+3;
    const double MHz=1.0E+6;
    const double mm=1.0E-3;
    const double cm=1.0E-2;

    // Definitions
    const int Ns=801;
    const int N_quadl=2;
    const double tol=1.0E-3;
    const int k_max=15;
    const double L=250.0*cm;
    const double h=3.0*cm;
    const double a=1.0*mm;
    const double freq_min=100.0*kHz;
    const double freq_max=100.0*MHz;
    const double lc=1.0/21.0;
    const cmplx V_1=1.0;
    const cmplx Z_1=50.0;
    const cmplx V_2=0.0;
    const cmplx Z_2=200.0;
    const double Z_c=(1.0/pi)*eta0*acosh(h/(2.0*a));

    // Program
    const cmplx j=cmplx(0.0, 1.0);
    double lambda;
    Range freq;    
    freq.logspace(freq_min, freq_max, Ns);
    File file;
    Engine engine;
    engine.quadl.set(N_quadl, tol, k_max);
    file.open("data/transmission_line/input_impedance/data1.dat", "w");
    cmplx Z_in, Z_in_ideal;
    for (int i=0; i<Ns; i++){
        std::cout << "Step " << i+1 << " out of " << Ns << std::endl;
        lambda = c0/freq(i);
        engine.shape.set_parameters(lambda, a, lc);
        engine.shape.create_transmission_line(L, h);
        engine.shape.add_port(1, V_1, Z_1);
        engine.shape.add_port(2, V_2, Z_2);
        engine.solve_delta_gap();
        Z_in = engine.find_Z_in(1);
        Z_in_ideal = Z_c*(Z_2+j*Z_c*tan(2.0*pi*L/lambda))/(Z_c+j*Z_2*tan(2.0*pi*L/lambda));
        file.write("%21.14E %21.14E %21.14E %21.14E %21.14E\n", freq(i), real(Z_in), imag(Z_in),
            real(Z_in_ideal), imag(Z_in_ideal));
        engine.deallocate();
    }
    file.close();

    T.unset();

}

void transmission_line_s_parameters(){

    Timer T;
    T.set();

    // Units
    const double MHz=1.0E+6;
    const double mm=1.0E-3;
    const double cm=1.0E-2;

    // Definitions
    const int Ns=801;
    const int N_quadl=2;
    const double tol=1.0E-3;
    const int k_max=15;
    const double L=250.0*cm;
    const double h=3.0*cm;
    const double a=1.0*mm;
    const double freq_min=0.1*MHz;
    const double freq_max=200.0*MHz;
    const double lc=1.0/21.0;
    const cmplx Z_0=50.0;
    const cmplx V_1=1.0;
    const cmplx Z_1=Z_0;
    const cmplx V_2=0.0;
    const cmplx Z_2=Z_0;
    const double Z_c=(1.0/pi)*eta0*acosh(h/(2.0*a));

    // Program
    const cmplx j=cmplx(0.0, 1.0);
    double lambda;
    Range freq;    
    freq.linspace(freq_min, freq_max, Ns);
    File file;
    Engine engine;
    engine.quadl.set(N_quadl, tol, k_max);
    file.open("data/transmission_line/s_parameters/data1.dat", "w");
    cmplx Z_in, Z_in_ideal;
    for (int i=0; i<Ns; i++){
        std::cout << "Step " << i+1 << " out of " << Ns << std::endl;
        lambda = c0/freq(i);
        engine.shape.set_parameters(lambda, a, lc);
        engine.shape.create_transmission_line(L, h);
        engine.shape.add_port(1, V_1, Z_1);
        engine.shape.add_port(2, V_2, Z_2);
        engine.solve_delta_gap();
        Z_in = engine.find_Z_in(1);
        int index=engine.shape.find_port_index(2);
        cmplx I_L=engine.I_n(index, 0);
        cmplx V_L=I_L*Z_0;        
        cmplx S_11=(Z_in-Z_0)/(Z_in+Z_0);
        cmplx S_21=2.0*V_L;
        cmplx beta_l=2.0*pi*L/lambda;
        Z_in_ideal = Z_c*(Z_2+j*Z_c*tan(beta_l))/(Z_c+j*Z_2*tan(beta_l));
        cmplx Gamma_l=(Z_0-Z_c)/(Z_0+Z_c);
        cmplx Gamma_s=(Z_c-Z_0)/(Z_c+Z_0);
        cmplx Gamma_0=Gamma_l*exp(-j*2.0*beta_l);
        cmplx S_11_ideal=(Gamma_s+Gamma_0)/(1.0+Gamma_s*Gamma_0);
        cmplx S_21_ideal=(1.0+Gamma_s)*(1.0+Gamma_l)*exp(-j*beta_l)/(1.0+Gamma_s*Gamma_0);;
        file.write("%21.14E %21.14E %21.14E %21.14E %21.14E\n", 
            freq(i)/MHz, 20*log10(abs(S_11)), 20*log10(abs(S_21)), 
            20*log10(abs(S_11_ideal)), 20*log10(abs(S_21_ideal)));
        engine.deallocate();
    }
    file.close();

    T.unset();

}

void yagi_uda_antenna_far_field(){

    Timer T;
    T.set();

    // Definitions
    const int Ns=801;
    const int N_quadl=2;
    const double tol=1.0E-3;
    const int k_max=15;
    const double lambda=1.0;
    const double a=0.003*lambda;
    const double lc=1.0/21.0;
    const cmplx V_in=1.0;
    const cmplx Z_in=50.0;
    const double theta_min=deg2rad(-180.0);
    const double theta_max=deg2rad(+180.0);
    const double phi_min=deg2rad(-180.0);
    const double phi_max=deg2rad(+180.0);

    Yag_Antenna yagi_antenna;
    yagi_antenna.N_reflectors = 1;
    yagi_antenna.N_directors = 13;
    double reflectors_lengths[1]={0.5*lambda};
    double reflectors_positions[1]={-0.25*lambda};
    double directos_lengths[13]={
        0.406*lambda, 0.406*lambda, 0.406*lambda, 0.406*lambda,
        0.406*lambda, 0.406*lambda, 0.406*lambda, 0.406*lambda,
        0.406*lambda, 0.406*lambda, 0.406*lambda, 0.406*lambda,
        0.406*lambda};
    double directos_positions[13];
    for (int i=0; i<13; i++){
        directos_positions[i] = (i+1)*0.34*lambda;
    }
    yagi_antenna.reflectors_legnths = reflectors_lengths;
    yagi_antenna.reflectors_positions = reflectors_positions;
    yagi_antenna.directors_legnths = directos_lengths;
    yagi_antenna.directors_positions = directos_positions;
    yagi_antenna.driver_length = 0.47*lambda;
    yagi_antenna.driver_position = 0.0*lambda;

    // Program
    Engine engine;
    engine.quadl.set(N_quadl, tol, k_max);
    engine.shape.set_parameters(lambda, a, lc);
    engine.shape.create_yagi_antenna(yagi_antenna);
    for (int i=0; i<13; i++){
        if (i==1){
            engine.shape.add_port(i+1, V_in, Z_in);
        }else{
            engine.shape.add_port(i+1, 0.0, 0.0);
        }
    }
    engine.shape.log("data/misc/shape_log.txt");
    engine.solve_delta_gap();
    Range theta, phi;    
    theta.linspace(theta_min, theta_max, Ns);
    phi.linspace(phi_min, phi_max, Ns);
    Far_Field far_field;
    File file;
    Matrix E_theta ,E_phi;
    E_theta.allocate(Ns, 1);
    E_phi.allocate(Ns, 1);
    // Elevation
    file.open("data/yagi_uda_antenna/far_field/data1.dat", "w");
    for (int i=0; i<Ns; i++){
        far_field = engine.get_far_field(theta(i), 0.0);
        E_theta(i, 0) = far_field.E_theta;
        E_phi(i, 0) = far_field.E_phi;
    }
    E_theta.normalize();
    E_phi.normalize();
    for (int i=0; i<Ns; i++){
        file.write("%21.14E %21.14E %21.14E\n", theta(i), 
            20.0*log10(abs(E_theta(i, 0))), 20.0*log10(abs(E_phi(i, 0))));
    }
    file.close();
    // Azimuth
    file.open("data/yagi_uda_antenna/far_field/data2.dat", "w");
    for (int i=0; i<Ns; i++){
        far_field = engine.get_far_field(pi/2.0, phi(i));
        E_theta(i, 0) = far_field.E_theta;
        E_phi(i, 0) = far_field.E_phi;
    }
    E_theta.normalize();
    E_phi.normalize();
    for (int i=0; i<Ns; i++){
        file.write("%21.14E %21.14E %21.14E\n", phi(i), 
            20.0*log10(abs(E_theta(i, 0))), 20.0*log10(abs(E_phi(i, 0))));
    }
    file.close();
    // Currents
    file.open("data/yagi_uda_antenna/far_field/data3.dat", "w");
    Matrix I_n;
    I_n.allocate(15, 1);
    for (int i=0; i<15; i++){
        int index=engine.shape.find_port_index(i+1);
        I_n(i, 0) = engine.I_n(index, 0);
    }
    I_n.normalize();
    for (int i=0; i<15; i++){
        file.write("%2d %21.14E\n", i+1, abs(I_n(i, 0)));
    }
    file.close();

    T.unset();

}

void vertical_dipole_RCS_1(){

    Timer T;
    T.set();

    // Units
    const double MHz=1.0E+6;

    // Definitions
    const int Ns=801;
    const int N_quadl=2;
    const double tol=1.0E-3;
    const int k_max=15;
    const double freq=600.0*MHz;
    const double lambda=c0/freq;
    const double L=1.0;
    const double a=L/200.0;
    const double lc=1.0/21.0;
    const double theta_i=deg2rad(45.0);
    const double phi_i=deg2rad(0.0);
    const double theta_min=deg2rad(0.0);
    const double theta_max=deg2rad(+180.0);
    const double phi=deg2rad(0.0);
    const cmplx E_TM=1.0;
    const cmplx E_TE=0.0;

    // Program
    Engine engine;
    engine.quadl.set(N_quadl, tol, k_max);
    engine.shape.set_parameters(lambda, a, lc);
    engine.shape.create_vertical_dipole(L, FALSE, 0);
    engine.solve_plane_wave(theta_i, phi_i, E_TM, E_TE);
    Range theta;    
    theta.linspace(theta_min, theta_max, Ns);
    RCS sigma;
    File file;
    // Elevation
    file.open("data/vertical_dipole/RCS/data1.dat", "w");
    for (int i=0; i<Ns; i++){
        sigma = engine.get_RCS(theta(i), phi);
        file.write("%21.14E %21.14E %21.14E\n", rad2deg(theta(i)), 
            10.0*log10(sigma.theta), 10.0*log10(sigma.phi));
    }
    file.close();

    T.unset();

}

void vertical_dipole_RCS_2(){

    Timer T;
    T.set();

    // Definitions
    const int Ns=801;
    const int N_quadl=2;
    const double tol=1.0E-3;
    const int k_max=15;
    const double lambda=1.0;
    const double lc=1.0/21.0;
    const double theta_i=deg2rad(30.0);
    const double phi_i=deg2rad(0.0);
    const double theta=deg2rad(30.0);
    const double phi=deg2rad(0.0);
    const double L_min=0.01*lambda;
    const double L_max=4.0*lambda;
    const cmplx E_TM=1.0;
    const cmplx E_TE=0.0;

    // Program
    Engine engine;
    engine.quadl.set(N_quadl, tol, k_max);
    Range L;    
    L.linspace(L_min, L_max, Ns);
    double a;
    RCS sigma;
    File file;
    // Elevation
    file.open("data/vertical_dipole/RCS/data2.dat", "w");
    for (int i=0; i<Ns; i++){
        std::cout << "Step " << i+1 << " out of " << Ns << std::endl;
        a = L(i)/200.0;
        engine.shape.set_parameters(lambda, a, lc);
        engine.shape.create_vertical_dipole(L(i), FALSE, 0);
        engine.solve_plane_wave(theta_i, phi_i, E_TM, E_TE);
        sigma = engine.get_RCS(theta, phi);
        file.write("%21.14E %21.14E %21.14E\n", L(i), 
            10.0*log10(sigma.theta), 10.0*log10(sigma.phi));
        engine.deallocate();
    }
    file.close();

    T.unset();

}

void vertical_dipole_RCS_3(){

    Timer T;
    T.set();

    // Definitions
    const int Ns=801;
    const int N_quadl=2;
    const double tol=1.0E-3;
    const int k_max=15;
    const double lambda=1.0;
    const double lc=1.0/21.0;
    const double theta_i=deg2rad(30.0);
    const double phi_i=deg2rad(0.0);
    const double theta=deg2rad(60.0);
    const double phi=deg2rad(0.0);
    const double L_min=0.01*lambda;
    const double L_max=4.0*lambda;
    const cmplx E_TM=1.0;
    const cmplx E_TE=0.0;

    // Program
    Engine engine;
    engine.quadl.set(N_quadl, tol, k_max);
    Range L;    
    L.linspace(L_min, L_max, Ns);
    double a;
    RCS sigma;
    File file;
    // Elevation
    file.open("data/vertical_dipole/RCS/data3.dat", "w");
    for (int i=0; i<Ns; i++){
        std::cout << "Step " << i+1 << " out of " << Ns << std::endl;
        a = L(i)/200.0;
        engine.shape.set_parameters(lambda, a, lc);
        engine.shape.create_vertical_dipole(L(i), FALSE, 0);
        engine.solve_plane_wave(theta_i, phi_i, E_TM, E_TE);
        sigma = engine.get_RCS(theta, phi);
        file.write("%21.14E %21.14E %21.14E\n", L(i), 
            10.0*log10(sigma.theta), 10.0*log10(sigma.phi));
        engine.deallocate();
    }
    file.close();

    T.unset();

}

void L_shape_RCS_1(){

    Timer T;
    T.set();

    // Units
    const double MHz=1.0E+6;

    // Definitions
    const int Ns=801;
    const int N_quadl=2;
    const double tol=1.0E-3;
    const int k_max=15;
    const double freq=600.0*MHz;
    const double lambda=c0/freq;
    const double L=1.0;
    const double a=L/200.0;
    const double lc=1.0/21.0;
    const double theta_i=deg2rad(60.0);
    const double phi_i=deg2rad(60.0);
    const double theta_min=deg2rad(0.0);
    const double theta_max=deg2rad(+180.0);
    const double phi=deg2rad(60.0);
    const cmplx E_TM=0.0;
    const cmplx E_TE=1.0;

    // Program
    Engine engine;
    engine.quadl.set(N_quadl, tol, k_max);
    engine.shape.set_parameters(lambda, a, lc);
    engine.shape.create_L_shape(L);
    engine.solve_plane_wave(theta_i, phi_i, E_TM, E_TE);
    Range theta;    
    theta.linspace(theta_min, theta_max, Ns);
    RCS sigma;
    File file;
    // Elevation
    file.open("data/L_shape/RCS/data1.dat", "w");
    for (int i=0; i<Ns; i++){
        sigma = engine.get_RCS(theta(i), phi);
        file.write("%21.14E %21.14E %21.14E\n", rad2deg(theta(i)), 
            10.0*log10(sigma.theta), 10.0*log10(sigma.phi));
    }
    file.close();

    T.unset();

}

void L_shape_RCS_2(){

    Timer T;
    T.set();

    // Units
    const double MHz=1.0E+6;

    // Definitions
    const int Ns=801;
    const int N_quadl=2;
    const double tol=1.0E-3;
    const int k_max=15;
    const double freq=900.0*MHz;
    const double lambda=c0/freq;
    const double L=1.0;
    const double a=L/200.0;
    const double lc=1.0/21.0;
    const double theta_i=deg2rad(30.0);
    const double phi_i=deg2rad(30.0);
    const double phi_min=deg2rad(-90.0);
    const double phi_max=deg2rad(+90.0);
    const double theta=deg2rad(60.0);
    const cmplx E_TM=1.0;
    const cmplx E_TE=0.0;

    // Program
    Engine engine;
    engine.quadl.set(N_quadl, tol, k_max);
    engine.shape.set_parameters(lambda, a, lc);
    engine.shape.create_L_shape(L);
    engine.solve_plane_wave(theta_i, phi_i, E_TM, E_TE);
    Range phi;    
    phi.linspace(phi_min, phi_max, Ns);
    RCS sigma;
    File file;
    // Elevation
    file.open("data/L_shape/RCS/data2.dat", "w");
    for (int i=0; i<Ns; i++){
        sigma = engine.get_RCS(theta, phi(i));
        file.write("%21.14E %21.14E %21.14E\n", rad2deg(phi(i)), 
            10.0*log10(sigma.theta), 10.0*log10(sigma.phi));
    }
    file.close();

    T.unset();

}

void monopole_s_parameters(){

    Timer T;
    T.set();

    // Units
    const double MHz=1.0E+6;
    const double GHz=1.0E+9;
    const double mm=1.0E-3;
    const double cm=1.0E-2;

    // Definitions
    const int Ns=801;
    const int N_quadl=2;
    const double tol=1.0E-3;
    const int k_max=15;
    const double L=10.5*cm;
    const double a=0.8*mm;
    const double freq_min=1.0*MHz;
    const double freq_max=10.0*GHz;
    const double lc=1.0/21.0;
    const cmplx Z_0=50.0;
    const cmplx V_1=1.0;
    const cmplx Z_1=Z_0;

    // Program
    double lambda;
    Range freq;    
    freq.linspace(freq_min, freq_max, Ns);
    File file;
    Engine engine;
    engine.quadl.set(N_quadl, tol, k_max);
    file.open("data/monopole/s_parameters/data1.dat", "w");
    cmplx Z_in, Z_in_ideal;
    for (int i=0; i<Ns; i++){
        std::cout << "Step " << i+1 << " out of " << Ns << std::endl;
        lambda = c0/freq(i);
        engine.shape.set_parameters(lambda, a, lc);
        engine.shape.create_vertical_dipole(2.0*L, TRUE, 1);
        engine.shape.add_port(1, 2.0*V_1, 2.0*Z_1);
        engine.solve_delta_gap();
        Z_in = engine.find_Z_in(1)/2.0;
        cmplx S_11=(Z_in-Z_0)/(Z_in+Z_0);
        file.write("%21.14E %21.14E\n", freq(i)/GHz, 20*log10(abs(S_11)));
        engine.deallocate();
    }
    file.close();

    T.unset();

}

void gull_antenna_far_field(){

    Timer T;
    T.set();

    // Units
    const double GHz=1.0E+9;
    const double mm=1.0E-3;

    // Definitions
    const int Ns=801;
    const int N_quadl=2;
    const double tol=1.0E-3;
    const int k_max=15;
    const double freq=3.0*GHz;
    const double lambda=c0/freq;
    const double h1=7.14*mm;
    const double h2=42.86*mm;
    const double h3=25.0*mm;
    const double alpha=deg2rad(50.0);
    const double a=0.5*mm;
    const double lc=1.0/21.0;
    const cmplx V_in=1.0;
    const cmplx Z_in=50.0;
    const double theta_min=deg2rad(-180.0);
    const double theta_max=deg2rad(+180.0);
    const double phi_min=deg2rad(-180.0);
    const double phi_max=deg2rad(+180.0);

    // Program
    Engine engine;
    engine.quadl.set(N_quadl, tol, k_max);
    engine.shape.set_parameters(lambda, a, lc);
    engine.shape.create_gull_antenna(h1, h2, h3, alpha);
    engine.shape.add_port(1, V_in, Z_in);
    engine.shape.log("data/misc/shape_log.txt");
    engine.solve_delta_gap();
    Range theta, phi;    
    theta.linspace(theta_min, theta_max, Ns);
    phi.linspace(phi_min, phi_max, Ns);
    Far_Field far_field;
    File file;
    Matrix E_theta ,E_phi;
    E_theta.allocate(Ns, 1);
    E_phi.allocate(Ns, 1);
    // Azimuth
    file.open("data/gull_antenna/far_field/data1.dat", "w");
    for (int i=0; i<Ns; i++){
        far_field = engine.get_far_field(theta(i), pi/2.0);
        E_theta(i, 0) = far_field.E_theta;
        E_phi(i, 0) = far_field.E_phi;
    }
    E_theta.normalize();
    E_phi.normalize();
    for (int i=0; i<Ns; i++){
        file.write("%21.14E %21.14E %21.14E\n", theta(i), 
            20.0*log10(abs(E_theta(i, 0))), 20.0*log10(abs(E_phi(i, 0))));
    }
    file.close();
    // Elevation
    file.open("data/gull_antenna/far_field/data2.dat", "w");
    for (int i=0; i<Ns; i++){
        far_field = engine.get_far_field(pi/2.0, phi(i));
        E_theta(i, 0) = far_field.E_theta;
        E_phi(i, 0) = far_field.E_phi;
    }
    E_theta.normalize();
    E_phi.normalize();
    for (int i=0; i<Ns; i++){
        file.write("%21.14E %21.14E %21.14E\n", phi(i), 
            20.0*log10(abs(E_theta(i, 0))), 20.0*log10(abs(E_phi(i, 0))));
    }
    file.close();

    T.unset();

}

void vertical_dipole_radiation_resistance(){

    Timer T;
    T.set();

    // Definitions
    const int Ns=801;
    const int N_quadl=2;
    const double tol=1.0E-3;
    const int k_max=15;
    const double lambda=1.0;
    const double a=1.0E-3*lambda;
    const double lc=1.0/21.0;
    const cmplx V_in=1.0;
    const cmplx Z_0=50.0;
    const double L_min=0.01*lambda;
    const double L_max=3.00*lambda;

    // Program
    Range L;    
    L.linspace(L_min, L_max, Ns);
    File file;
    Engine engine;
    engine.quadl.set(N_quadl, tol, k_max);
    file.open("data/vertical_dipole/radiation_resistance/data1.dat", "w");
    cmplx Z_in;
    double R_rad;
    for (int i=0; i<Ns; i++){
        std::cout << "Step " << i+1 << " out of " << Ns << std::endl;
        engine.shape.set_parameters(lambda, a, lc);
        engine.shape.create_vertical_dipole(L(i), TRUE, 1);
        engine.shape.add_port(1, V_in, Z_0);
        engine.solve_delta_gap();
        Z_in = engine.find_Z_in(1);
        R_rad = engine.radiation_resistance(1);
        file.write("%21.14E %21.14E %21.14E\n", L(i)/lambda, real(Z_in), R_rad);
        engine.deallocate();
    }
    file.close();

    T.unset();

}

void vertical_dipole_directive_gain(){

    Timer T;
    T.set();

    // Definitions
    const int Ns=801;
    const int N_quadl=2;
    const double tol=1.0E-3;
    const int k_max=15;
    const double lambda=1.0;
    const double L=0.47*lambda;
    const double a=1.0E-3*lambda;
    const double lc=1.0/21.0;
    const cmplx V_in=1.0;
    const cmplx Z_0=50.0;
    const double theta_min=deg2rad(-180.0);
    const double theta_max=deg2rad(+180.0);
    const double phi_min=deg2rad(-180.0);
    const double phi_max=deg2rad(+180.0);

    // Program
    Engine engine;
    engine.quadl.set(N_quadl, tol, k_max);
    engine.shape.set_parameters(lambda, a, lc);
    engine.shape.create_vertical_dipole(L, TRUE, 1);
    engine.shape.add_port(1, V_in, Z_0);
    engine.solve_delta_gap();
    Range theta, phi;    
    theta.linspace(theta_min, theta_max, Ns);
    phi.linspace(phi_min, phi_max, Ns);
    Far_Field far_field;
    File file;
    Matrix directive_gain;
    directive_gain.allocate(Ns, 1);
    // Elevation
    file.open("data/vertical_dipole/directive_gain/data1.dat", "w");
    for (int i=0; i<Ns; i++){
        directive_gain(i, 0) = engine.directive_gain(theta(i), 0.0);
        file.write("%21.14E %21.14E\n", theta(i), 
            10.0*log10(abs(directive_gain(i, 0))));
        progress_bar(i, Ns, "computing directive gain...");
    }
    file.close();
    // Azimuth
    file.open("data/vertical_dipole/directive_gain/data2.dat", "w");
    for (int i=0; i<Ns; i++){
        directive_gain(i, 0) = engine.directive_gain(pi/2.0, phi(i));
        file.write("%21.14E %21.14E\n", theta(i), 
            10.0*log10(abs(directive_gain(i, 0))));
        progress_bar(i, Ns, "computing directive gain...");
    }
    file.close();

    T.unset();

}

void transmission_line_near_field(){

    Timer T;
    T.set();

    // Units
    const double GHz=1.0E+9;
    const double mm=1.0E-3;
    const double cm=1.0E-2;
    const double mA=1.0E-3;

    double factor=1.0;

    // Definitions
    const int Ns=51;
    const int N_quadl=2;
    const double tol=1.0E-3;
    const int k_max=15;
    const double freq=3.0*GHz*factor;
    const double lambda=c0/freq;
    const double L=10*cm/factor;
    const double h=2.0*5*cm/factor;
    const double a=0.8*mm/factor;
    const double lc=1.0/21.0;
    const cmplx V_1=1.0;
    const cmplx Z_1=50.0;
    const cmplx V_2=0.0;
    const cmplx Z_2=200.0;
    const double x=+10.0*cm/factor;
    const double y_min=-10.0*cm/factor;
    const double y_max=+10.0*cm/factor;
    const double z=+2.0*cm/factor;

    // Program
    Engine engine;
    engine.quadl.set(N_quadl, tol, k_max);
    engine.shape.set_parameters(lambda, a, lc);
    engine.shape.create_transmission_line(L, h);
    engine.shape.add_port(1, 2.0*V_1, 2.0*Z_1);
    engine.shape.add_port(2, 2.0*V_2, 2.0*Z_2);
    engine.solve_delta_gap();
    Range y;    
    y.linspace(y_min, y_max, Ns);
    Near_Field near_field;
    File file_E, file_H;
    // E field
    file_E.open("data/transmission_line/near_field/data1.dat", "w");
    file_H.open("data/transmission_line/near_field/data2.dat", "w");
    double E_total, H_total;
    for (int i=0; i<Ns; i++){
        near_field = engine.get_near_field(x, y(i), z);
        E_total = sqrt(
            pow(abs(near_field.E_x), 2.0)+
            pow(abs(near_field.E_y), 2.0)+
            pow(abs(near_field.E_z), 2.0)
        );
        H_total = sqrt(
            pow(abs(near_field.H_x), 2.0)+
            pow(abs(near_field.H_y), 2.0)+
            pow(abs(near_field.H_z), 2.0)
        );
        file_E.write("%21.14E %21.14E %21.14E %21.14E %21.14E\n", y(i)/cm, 
            abs(near_field.E_x), abs(near_field.E_y), abs(near_field.E_z), E_total);
        file_H.write("%21.14E %21.14E %21.14E %21.14E %21.14E\n", y(i)/cm, 
            abs(near_field.H_x)/mA, abs(near_field.H_y)/mA, abs(near_field.H_z)/mA, H_total);
        progress_bar(i, Ns, "computing near field...");
    }
    file_E.close();
    file_H.close();
    disp(engine.find_Z_in(1)/2.0);

    T.unset();

}

void vertical_dipole_near_field_2D(){

    Timer T;
    T.set();

    // Definitions
    const int Ns_x=801;
    const int Ns_z=801;
    const int N_quadl=2;
    const double tol=1.0E-3;
    const int k_max=15;
    const double lambda=3.0;
    const double L=0.47*lambda;
    const double a=5.0E-3*lambda;
    const double lc=1.0/21.0;
    const cmplx V_in=1.0;
    const cmplx Z_0=50.0;
    const double y=0.0*lambda;
    const double x_min=-4.0*lambda;
    const double x_max=+4.0*lambda;
    const double z_min=-4.0*lambda;
    const double z_max=+4.0*lambda;

    // Program
    Range x, z;    
    x.linspace(x_min, x_max, Ns_x);
    z.linspace(z_min, z_max, Ns_z);
    File file;
    Engine engine;
    engine.quadl.set(N_quadl, tol, k_max);
    engine.shape.set_parameters(lambda, a, lc);
    engine.shape.create_vertical_dipole(L, TRUE, 1);
    engine.shape.add_port(1, V_in, Z_0);
    engine.solve_delta_gap();
    engine.I_n.save("data/misc/I_n.dat");
    engine.shape.log("data/misc/shape_log.txt");
    file.open("data/vertical_dipole/near_field_2D/data1.dat", "w");
    int counter=0;
    Near_Field near_field;
    for (int i=0; i<Ns_x; i++){
        for (int j=0; j<Ns_z; j++){
            near_field = engine.get_near_field(x(i), y, z(j));
            double E=sqrt(pow(real(near_field.E_x), 2.0)+
                          pow(real(near_field.E_y), 2.0)+
                          pow(real(near_field.E_z), 2.0));
            file.write("%21.14E %21.14E %21.14E\n", x(i)/lambda, z(j)/lambda, 20.0*log10(abs(E)));
            progress_bar(counter++, Ns_x*Ns_z, "computing near fields...");
        }
        file.write("\n");
    }
    file.close();

    T.unset();

}

void vertical_dipole_near_field_vs_far_field(){

    Timer T;
    T.set();

    // Definitions
    const int Ns=801;
    const int N_quadl=2;
    const double tol=1.0E-3;
    const int k_max=15;
    const double lambda=1.0;
    const double L=0.47*lambda;
    const double a=5.0E-3*lambda;
    const double lc=1.0/21.0;
    const cmplx V_in=1.0;
    const cmplx Z_0=50.0;
    const double x_min=0.01*lambda;
    const double x_max=100.0*lambda;

    // Program
    Range x;    
    x.logspace(x_min, x_max, Ns);
    File file_E, file_H;
    Engine engine;
    engine.quadl.set(N_quadl, tol, k_max);
    engine.shape.set_parameters(lambda, a, lc);
    engine.shape.create_vertical_dipole(L, TRUE, 1);
    engine.shape.add_port(1, V_in, Z_0);
    engine.solve_delta_gap();
    file_E.open("data/vertical_dipole/near_vs_far_field/data1.dat", "w");
    file_H.open("data/vertical_dipole/near_vs_far_field/data2.dat", "w");
    const cmplx j=cmplx(0.0, 1.0);
    const double k=2.0*pi;
    cmplx factor;
    Near_Field near_field;
    Far_Field far_field;
    double E_near, E_far;
    double H_near, H_far;
    for (int i=0; i<Ns; i++){
        factor = exp(-j*k*x(i))/(4.0*pi*x(i));
        near_field = engine.get_near_field(x(i), 0.0, 0.0);
        far_field = engine.get_far_field(pi/2.0, 0.0);
        E_near = sqrt(pow(abs(near_field.E_x), 2.0)+
                      pow(abs(near_field.E_y), 2.0)+
                      pow(abs(near_field.E_z), 2.0));
        E_far = abs(factor)*sqrt(pow(abs(far_field.E_theta), 2.0)+
                                 pow(abs(far_field.E_phi), 2.0));  
        H_near = sqrt(pow(abs(near_field.H_x), 2.0)+
                      pow(abs(near_field.H_y), 2.0)+
                      pow(abs(near_field.H_z), 2.0));
        H_far = abs(factor)*sqrt(pow(abs(far_field.H_theta), 2.0)+
                                 pow(abs(far_field.H_phi), 2.0)); 
        file_E.write("%21.14E %21.14E %21.14E\n", x(i), 20.0*log10(E_near), 20.0*log10(E_far));
        file_H.write("%21.14E %21.14E %21.14E\n", x(i), 20.0*log10(H_near), 20.0*log10(H_far));
        progress_bar(i, Ns, "computing fields...");
    }
    file_E.close();
    file_H.close();

    T.unset();

}

void transmission_line_far_field(){

    Timer T;
    T.set();

    // Units
    const double GHz=1.0E+9;
    const double mm=1.0E-3;
    const double cm=1.0E-2;

    // Definitions
    const int Ns=801;
    const int N_quadl=2;
    const double tol=1.0E-3;
    const int k_max=15;
    const double freq=1.0*GHz;
    const double lambda=c0/freq;
    const double L=150.0*cm;
    const double h=50.0*cm;
    const double a=0.8*mm;
    const double lc=1.0/21.0;
    const cmplx V_1=1.0;
    const cmplx Z_1=50.0;
    const cmplx V_2=0.0;
    const cmplx Z_2=0.0;
    const double theta_min=deg2rad(-180.0);
    const double theta_max=deg2rad(+180.0);

    // Program
    Range theta;    
    theta.linspace(theta_min, theta_max, Ns);
    File file;
    Engine engine;
    engine.quadl.set(N_quadl, tol, k_max);
    engine.shape.set_parameters(lambda, a, lc);
    engine.shape.create_transmission_line(L, h);
    engine.shape.add_port(1, 2.0*V_1, 2.0*Z_1);
    engine.shape.add_port(2, 2.0*V_2, 2.0*Z_2);
    engine.solve_delta_gap();
    file.open("data/transmission_line/far_field/data1.dat", "w");
    Far_Field far_field;
    Matrix E_theta;
    E_theta.allocate(Ns, 1);
    for (int i=0; i<Ns; i++){
        far_field = engine.get_far_field(theta(i), 0.0);
        E_theta(i, 0) = far_field.E_theta;
        progress_bar(i, Ns, "computing far fields...");
    }
    E_theta.normalize();
    for (int i=0; i<Ns; i++){
        disp(theta(i));
        file.write("%21.14E %21.14E\n", theta(i), 20.0*log10(abs(E_theta(i, 0))));
    }
    file.close();

    T.unset();

}

void RF_coil_s_parameters(){

    Timer T;
    T.set();

    // Units
    const double MHz=1.0E+6;
    const double cm=1.0E-2;
    const double mm=1.0E-3;

    // Definitions
    const int Ns=401;
    const int N_quadl=2;
    const double tol=1.0E-3;
    const int k_max=15;
    const double R=50.0*cm;
    const double H=100.0*cm;
    const double a=2.0*mm/(pi);
    const int N_legs=8;
    const double freq_min=10.0*MHz;
    const double freq_max=100.0*MHz;
    const double lc=1.0/21.0;
    const cmplx Z_0=50.0;
    const int N_ports=N_legs;

    // Program
    double lambda;
    Range freq;    
    freq.linspace(freq_min, freq_max, Ns);
    File file;
    Engine engine;
    engine.quadl.set(N_quadl, tol, k_max);
    file.open("data/RF_coil/s_parameters/data.csv", "w");
    cmplx Z_in, Z_in_ideal;
    Matrix S;
    S.allocate(N_ports, N_ports);
    for (int i=0; i<Ns; i++){
        std::cout << "Step " << i+1 << " out of " << Ns << std::endl;
        lambda = c0/freq(i);
        engine.shape.set_parameters(lambda, a, lc);
        engine.shape.create_RF_coil(R, H, N_legs);
        engine.get_S_paramters(Z_0, S);
        file.write("%21.14E, ", freq(i));
        for (int m=0; m<N_legs; m++){
            for (int n=0; n<N_legs; n++){
                if (m==N_legs-1 && n==N_legs-1){
                    file.write("%21.14E, ", real(S(m, n)));
                    file.write("%21.14E", imag(S(m, n)));
                }else{
                    file.write("%21.14E, ", real(S(m, n)));
                    file.write("%21.14E, ", imag(S(m, n)));
                }
            }
        }
        file.write("\n");
        engine.deallocate();
    }
    file.close();

    T.unset();

}