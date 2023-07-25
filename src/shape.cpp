//
#include "shape.hpp"

// Shape class
Shape::Shape(){}

Shape::~Shape(){
    Shape::deallocate();
}

void Shape::deallocate(){
    if (this->is_allocated){
        free(this->bases_list);
    }
    this->is_allocated = FALSE;
    this->is_set = FALSE;
    this->N_bases = 0;
    this->lambda = 0.0;
    this->a = 0.0;
    this->lc = 0.0;
    this->N_ports = 0;
}

void Shape::set_parameters(const double lambda, const double a, const double lc){
    check_error(this->is_set, "shape parameters are already set!");
    check_error(lambda<=0.0, "invalid value for lambda!");
    check_error(a<=0.0, "invalid value for wire radius!");
    check_error(lambda<=0.0, "invalid lc value!");
    this->lambda = lambda;
    this->a = a/this->lambda;
    this->lc = lc;
    this->is_set = TRUE;
}

int Shape::get_N_ports(){
    this->N_ports = 0;
    for (int i=0; i<this->N_bases; i++){
        if (this->bases_list[i].port_number>0){
            this->N_ports++;
        }
    }
    return this->N_ports;
}

void Shape::allocate(){
    File file;
    file.open("mesh/shape_info.txt", "r");
    file.read("%d", &this->N_bases);
    check_error(this->N_bases<1, "invalid mesh data!");
    file.close();
    file.open("mesh/bases_info.txt", "r");
    this->bases_list = (Basis*)calloc(this->N_bases, sizeof(Basis));
	assert(this->bases_list!=NULL);
    for (int i=0; i<this->N_bases; i++){
        Basis basis;
        file.read("%lf %lf %lf", &basis.r_m.x, &basis.r_m.y, &basis.r_m.z);
        file.read("%lf %lf %lf", &basis.r_n.x, &basis.r_n.y, &basis.r_n.z);
        file.read("%lf %lf %lf", &basis.r_p.x, &basis.r_p.y, &basis.r_p.z);
        file.read("%d", &basis.port_number);
        basis.index = i;
        basis.get_parameters();
        this->bases_list[i] = basis;
    }
    file.close();
    this->N_ports = 0;
    this->is_allocated = TRUE;
}

void Shape::get_parameters(int &N_bases, double &lambda, double &a){
    N_bases = this->N_bases;
    lambda = this->lambda;
    a = this->a;
}

void Shape::add_port(const int port_number, const cmplx V, const cmplx Z){
    check_error(!this->is_allocated, "shape is not allocated yet!");
    int index=Shape::find_port_index(port_number);
    check_error(index<0, "invalid port number!");
    this->bases_list[index].V = V;
    this->bases_list[index].Z = Z;
    this->N_ports++;
}

int Shape::find_port_index(const int port_number){
    check_error(port_number<1, "invalid port number!");
    for (int i=0; i<this->N_bases; i++){
        if (this->bases_list[i].port_number==port_number){
            return this->bases_list[i].index;
        }
    }
    check_error(TRUE, "invalid port number!");
    return -1;
}

Basis* Shape::get_basis(const int index){
    check_error(index<0, "invalid basis index!");
    check_error(index>=this->N_bases, "invalid basis index!");
    return &(this->bases_list[index]);
}

void Shape::reset_ports(){
    for (int i=0; i<this->N_bases; i++){
        this->bases_list[i].V = 0.0;
        this->bases_list[i].Z = 0.0;
    }
}

void Shape::log(const char *filename){
    File file;
    file.open(filename, "w");
    file.write("wavelength: %21.14E\n", this->lambda);
    file.write("wire radius: %21.14E\n", this->a);
    file.write("lc: %21.14E\n", this->lc);
    file.write("number of basis functions: %d\n\n", this->N_bases);
    for (int i=0; i<this->N_bases; i++){
        file.write("basis %d:\n", this->bases_list[i].index);
        file.write("port number: %d\n", this->bases_list[i].port_number);
        file.write("voltage: (%21.14E, %21.14E)\n", 
            real(this->bases_list[i].V),
            imag(this->bases_list[i].V));
        file.write("impedance: (%21.14E, %21.14E)\n", 
            real(this->bases_list[i].Z),
            imag(this->bases_list[i].Z));
        file.write("r_m.x: (%21.14E, %21.14E)\n", 
            real(this->bases_list[i].r_m.x),
            imag(this->bases_list[i].r_m.x));
        file.write("r_m.y: (%21.14E, %21.14E)\n", 
            real(this->bases_list[i].r_m.y),
            imag(this->bases_list[i].r_m.y));
        file.write("r_m.z: (%21.14E, %21.14E)\n", 
            real(this->bases_list[i].r_m.z),
            imag(this->bases_list[i].r_m.z));
        file.write("r_n.x: (%21.14E, %21.14E)\n", 
            real(this->bases_list[i].r_n.x),
            imag(this->bases_list[i].r_n.x));
        file.write("r_n.y: (%21.14E, %21.14E)\n", 
            real(this->bases_list[i].r_n.y),
            imag(this->bases_list[i].r_n.y));
        file.write("r_n.z: (%21.14E, %21.14E)\n", 
            real(this->bases_list[i].r_n.z),
            imag(this->bases_list[i].r_n.z));
        file.write("r_p.x: (%21.14E, %21.14E)\n", 
            real(this->bases_list[i].r_p.x),
            imag(this->bases_list[i].r_p.x));
        file.write("r_p.y: (%21.14E, %21.14E)\n", 
            real(this->bases_list[i].r_p.y),
            imag(this->bases_list[i].r_p.y));
        file.write("r_p.z: (%21.14E, %21.14E)\n", 
            real(this->bases_list[i].r_p.z),
            imag(this->bases_list[i].r_p.z));
        file.write("L_m.x: (%21.14E, %21.14E)\n", 
            real(this->bases_list[i].L_m.x),
            imag(this->bases_list[i].L_m.x));
        file.write("L_m.y: (%21.14E, %21.14E)\n", 
            real(this->bases_list[i].L_m.y),
            imag(this->bases_list[i].L_m.y));
        file.write("L_m.z: (%21.14E, %21.14E)\n", 
            real(this->bases_list[i].L_m.z),
            imag(this->bases_list[i].L_m.z));
        file.write("L_p.x: (%21.14E, %21.14E)\n", 
            real(this->bases_list[i].L_p.x),
            imag(this->bases_list[i].L_p.x));
        file.write("L_p.y: (%21.14E, %21.14E)\n", 
            real(this->bases_list[i].L_p.y),
            imag(this->bases_list[i].L_p.y));
        file.write("L_p.z: (%21.14E, %21.14E)\n", 
            real(this->bases_list[i].L_p.z),
            imag(this->bases_list[i].L_p.z));
        file.write("l_m: %21.14E\n", this->bases_list[i].l_m);
        file.write("l_p: %21.14E\n\n", this->bases_list[i].l_p);
    }
    file.close();
}

// Shapes library
void Shape::create_vertical_dipole(const double L, const int is_port, const int port_number){
    // 
    check_error(this->is_allocated, "shape is already allocated!");
    check_error(!this->is_set, "shape parameters are not set yet!");
    File file;
    file.open("mesh/shape.py", "w");
    file.write("import numpy as np\n");
    file.write("import mesh\n");
    file.write("lambda_ = %21.14E\n", this->lambda);
    file.write("lc = %21.14E\n", this->lc);
    file.write("lines_list = []\n");
    file.write("segments_list = []\n");
    // Begin shape definition
    check_error(L<=0.0, "invalid dipole length!");
    file.write("L = %21.14E\n", L);
    file.write("v1 = mesh.Vertex(0.0, 0.0, -L/2.0)\n");
    file.write("v2 = mesh.Vertex(0.0, 0.0, +L/2.0)\n");
    if (is_port){
        file.write("lines_list.append(mesh.Line(v1, v2, True, %d))\n", port_number);
    }else{
        file.write("lines_list.append(mesh.Line(v1, v2, False, 0))\n", 0);
    }
    // End of shape definition
    file.write("shape = mesh.Shape(lambda_, lc)\n");
    file.write("for line in lines_list:\n");
    file.write("\tshape.add_line(line)\n");
    file.write("\tcontinue\n");
    file.write("for segment in segments_list:\n");
    file.write("\tshape.add_segment(segment)\n");
    file.write("\tcontinue\n");
    file.write("shape.mesh()\n");
    file.close();
#ifdef _WIN64
    assert(!system("python mesh/shape.py"));
#endif
#ifdef __linux__
    assert(!system("python3 mesh/shape.py"));
#endif
    Shape::allocate();
}

void Shape::create_circular_loop(const double r, const int is_port, const int port_number){
    // 
    check_error(this->is_allocated, "shape is already allocated!");
    check_error(!this->is_set, "shape parameters are not set yet!");
    File file;
    file.open("mesh/shape.py", "w");
    file.write("import numpy as np\n");
    file.write("import mesh\n");
    file.write("lambda_ = %21.14E\n", this->lambda);
    file.write("lc = %21.14E\n", this->lc);
    file.write("lines_list = []\n");
    file.write("segments_list = []\n");
    // Begin shape definition
    check_error(r<=0.0, "invalid loop radius!");
    file.write("r = %21.14E\n", r);
    file.write("N = int(np.ceil(2.0*np.pi*(r/lambda_)/lc))\n");
    file.write("if N<6:\n");
    file.write("\tN = 6\n");
    file.write("\tpass\n");
    file.write("d_phi = 2.0*np.pi/N\n");
    file.write("for i in range(N):\n");
    file.write("\tphi = (i+0)*d_phi\n");
    file.write("\tv1 = mesh.Vertex(r*np.cos(phi), r*np.sin(phi), 0.0)\n");
    file.write("\tphi = (i+1)*d_phi\n");
    file.write("\tv2 = mesh.Vertex(r*np.cos(phi), r*np.sin(phi), 0.0)\n");
    file.write("\tif i==0 or i==N-1:\n");
    if (is_port){
        file.write("\t\tsegment = mesh.Segment(v1, v2, True, %d)\n", port_number);
    }else{
        file.write("\t\tsegment = mesh.Segment(v1, v2, False, 0)\n");    
    }
    file.write("\t\tpass\n");
    file.write("\telse:\n");
    file.write("\t\tsegment = mesh.Segment(v1, v2, False, 0)\n");
    file.write("\t\tpass\n");
    file.write("\tsegments_list.append(segment)\n");
    file.write("\tcontinue\n");
    // End of shape definition
    file.write("shape = mesh.Shape(lambda_, lc)\n");
    file.write("for line in lines_list:\n");
    file.write("\tshape.add_line(line)\n");
    file.write("\tcontinue\n");
    file.write("for segment in segments_list:\n");
    file.write("\tshape.add_segment(segment)\n");
    file.write("\tcontinue\n");
    file.write("shape.mesh()\n");
    file.close();
#ifdef _WIN64
    assert(!system("python mesh/shape.py"));
#endif
#ifdef __linux__
    assert(!system("python3 mesh/shape.py"));
#endif
    Shape::allocate();
}

void Shape::create_transmission_line(const double L, const double h){
    // 
    check_error(this->is_allocated, "shape is already allocated!");
    check_error(!this->is_set, "shape parameters are not set yet!");
    File file;
    file.open("mesh/shape.py", "w");
    file.write("import numpy as np\n");
    file.write("import mesh\n");
    file.write("lambda_ = %21.14E\n", this->lambda);
    file.write("lc = %21.14E\n", this->lc);
    file.write("lines_list = []\n");
    file.write("segments_list = []\n");
    // Begin shape definition
    check_error(L<=0.0, "invalid transmission line length!");
    check_error(h<=0.0, "invalid transmission line height!");
    file.write("L = %21.14E\n", L);
    file.write("h = %21.14E\n", h);
    file.write("v1 = mesh.Vertex(-L/2.0, 0.0, -h/2.0)\n");
    file.write("v2 = mesh.Vertex(-L/2.0, 0.0, +h/2.0)\n");
    file.write("v3 = mesh.Vertex(+L/2.0, 0.0, +h/2.0)\n");
    file.write("v4 = mesh.Vertex(+L/2.0, 0.0, -h/2.0)\n");
    file.write("lines_list.append(mesh.Line(v1, v2, True, 1))\n");
    file.write("lines_list.append(mesh.Line(v2, v3, False, 0))\n");
    file.write("lines_list.append(mesh.Line(v4, v3, True, 2))\n");
    file.write("lines_list.append(mesh.Line(v4, v1, False, 0))\n");
    // End of shape definition
    file.write("shape = mesh.Shape(lambda_, lc)\n");
    file.write("for line in lines_list:\n");
    file.write("\tshape.add_line(line)\n");
    file.write("\tcontinue\n");
    file.write("for segment in segments_list:\n");
    file.write("\tshape.add_segment(segment)\n");
    file.write("\tcontinue\n");
    file.write("shape.mesh()\n");
    file.close();
#ifdef _WIN64
    assert(!system("python mesh/shape.py"));
#endif
#ifdef __linux__
    assert(!system("python3 mesh/shape.py"));
#endif
    Shape::allocate();
}

void Shape::create_yagi_antenna(const Yag_Antenna yagi_antnna){
    // 
    check_error(this->is_allocated, "shape is already allocated!");
    check_error(!this->is_set, "shape parameters are not set yet!");
    File file;
    file.open("mesh/shape.py", "w");
    file.write("import numpy as np\n");
    file.write("import mesh\n");
    file.write("lambda_ = %21.14E\n", this->lambda);
    file.write("lc = %21.14E\n", this->lc);
    file.write("lines_list = []\n");
    file.write("segments_list = []\n");
    // Begin shape definition
    check_error(yagi_antnna.directors_legnths==NULL, "yagi directors info are missing!");
    check_error(yagi_antnna.directors_positions==NULL, "yagi directors info are missing!");
    check_error(yagi_antnna.reflectors_legnths==NULL, "yagi reflectors info are missing!");
    check_error(yagi_antnna.reflectors_positions==NULL, "yagi reflectors info are missing!");
    double L;
    double S;
    int port_number=0;
    // Reflectors
    int N_reflectors=yagi_antnna.N_reflectors;
    check_error(N_reflectors<1, "invalid yagi antenna!");
    for (int i=0; i<N_reflectors; i++){
        port_number++;
        S = yagi_antnna.reflectors_positions[i];
        L = yagi_antnna.reflectors_legnths[i];
        check_error(L<=0.0, "invalid yagi antenna!");
        file.write("v1 = mesh.Vertex(%21.14E, 0.0, %21.14E)\n", S, -L/2.0);
        file.write("v2 = mesh.Vertex(%21.14E, 0.0, %21.14E)\n", S, +L/2.0);
        file.write("lines_list.append(mesh.Line(v1, v2, True, %d))\n", port_number);
    }
    // Driver
    port_number++;
    S = yagi_antnna.driver_position;
    L = yagi_antnna.driver_length;
    file.write("v1 = mesh.Vertex(%21.14E, 0.0, %21.14E)\n", S, -L/2.0);
    file.write("v2 = mesh.Vertex(%21.14E, 0.0, %21.14E)\n", S, +L/2.0);
    file.write("lines_list.append(mesh.Line(v1, v2, True, %d))\n", port_number);
    // Directors
    int N_directors=yagi_antnna.N_directors;
    check_error(N_directors<1, "invalid yagi antenna!");
    for (int i=0; i<N_directors; i++){
        port_number++;
        S = yagi_antnna.directors_positions[i];
        L = yagi_antnna.directors_legnths[i];
        check_error(L<=0.0, "invalid yagi antenna!");
        file.write("v1 = mesh.Vertex(%21.14E, 0.0, %21.14E)\n", S, -L/2.0);
        file.write("v2 = mesh.Vertex(%21.14E, 0.0, %21.14E)\n", S, +L/2.0);
        file.write("lines_list.append(mesh.Line(v1, v2, True, %d))\n", port_number);
    }
    // End of shape definition
    file.write("shape = mesh.Shape(lambda_, lc)\n");
    file.write("for line in lines_list:\n");
    file.write("\tshape.add_line(line)\n");
    file.write("\tcontinue\n");
    file.write("for segment in segments_list:\n");
    file.write("\tshape.add_segment(segment)\n");
    file.write("\tcontinue\n");
    file.write("shape.mesh()\n");
    file.close();
#ifdef _WIN64
    assert(!system("python mesh/shape.py"));
#endif
#ifdef __linux__
    assert(!system("python3 mesh/shape.py"));
#endif
    Shape::allocate();
}

void Shape::create_L_shape(const double L){
    // 
    check_error(this->is_allocated, "shape is already allocated!");
    check_error(!this->is_set, "shape parameters are not set yet!");
    File file;
    file.open("mesh/shape.py", "w");
    file.write("import numpy as np\n");
    file.write("import mesh\n");
    file.write("lambda_ = %21.14E\n", this->lambda);
    file.write("lc = %21.14E\n", this->lc);
    file.write("lines_list = []\n");
    file.write("segments_list = []\n");
    // Begin shape definition
    check_error(L<=0.0, "invalid length!");
    file.write("L = %21.14E\n", L);
    file.write("v1 = mesh.Vertex(0.0, 0.0, 0.0)\n");
    file.write("v2 = mesh.Vertex(0.0, 0.0, 0.3*L)\n");
    file.write("v3 = mesh.Vertex(0.0, 0.7*L, 0.0)\n");
    file.write("lines_list.append(mesh.Line(v1, v2, False, 0))\n");
    file.write("lines_list.append(mesh.Line(v1, v3, False, 0))\n");
    // End of shape definition
    file.write("shape = mesh.Shape(lambda_, lc)\n");
    file.write("for line in lines_list:\n");
    file.write("\tshape.add_line(line)\n");
    file.write("\tcontinue\n");
    file.write("for segment in segments_list:\n");
    file.write("\tshape.add_segment(segment)\n");
    file.write("\tcontinue\n");
    file.write("shape.mesh()\n");
    file.close();
#ifdef _WIN64
    assert(!system("python mesh/shape.py"));
#endif
#ifdef __linux__
    assert(!system("python3 mesh/shape.py"));
#endif
    Shape::allocate();
}

void Shape::create_gull_antenna(const double h1, const double h2, 
        const double h3, const double alpha){
    // 
    check_error(this->is_allocated, "shape is already allocated!");
    check_error(!this->is_set, "shape parameters are not set yet!");
    File file;
    file.open("mesh/shape.py", "w");
    file.write("import numpy as np\n");
    file.write("import mesh\n");
    file.write("lambda_ = %21.14E\n", this->lambda);
    file.write("lc = %21.14E\n", this->lc);
    file.write("lines_list = []\n");
    file.write("segments_list = []\n");
    // Begin shape definition
    check_error(h1<=0.0||h2<=0.0||h3<=0.0, "invalid length!");
    file.write("h1 = %21.14E\n", h1);
    file.write("h2 = %21.14E\n", h2);
    file.write("h3 = %21.14E\n", h3);
    file.write("alpha = %21.14E\n", alpha);
    file.write("v1 = mesh.Vertex(-h1, 0.0, 0.0)\n");
    file.write("v2 = mesh.Vertex(+h1, 0.0, 0.0)\n");
    file.write("lines_list.append(mesh.Line(v1, v2, True, 1))\n");
    file.write("v3 = mesh.Vertex(-h1-h2*np.cos(alpha), +h2*np.sin(alpha), 0.0)\n");
    file.write("v4 = mesh.Vertex(+h1+h2*np.cos(alpha), +h2*np.sin(alpha), 0.0)\n");
    file.write("v5 = mesh.Vertex(-h1-h2*np.cos(alpha)-h3, +h2*np.sin(alpha), 0.0)\n");
    file.write("v6 = mesh.Vertex(+h1+h2*np.cos(alpha)+h3, +h2*np.sin(alpha), 0.0)\n");
    file.write("lines_list.append(mesh.Line(v1, v3, False, 0))\n");
    file.write("lines_list.append(mesh.Line(v3, v5, False, 0))\n");
    file.write("lines_list.append(mesh.Line(v2, v4, False, 0))\n");
    file.write("lines_list.append(mesh.Line(v4, v6, False, 0))\n");
    // End of shape definition
    file.write("shape = mesh.Shape(lambda_, lc)\n");
    file.write("for line in lines_list:\n");
    file.write("\tshape.add_line(line)\n");
    file.write("\tcontinue\n");
    file.write("for segment in segments_list:\n");
    file.write("\tshape.add_segment(segment)\n");
    file.write("\tcontinue\n");
    file.write("shape.mesh()\n");
    file.close();
#ifdef _WIN64
    assert(!system("python mesh/shape.py"));
#endif
#ifdef __linux__
    assert(!system("python3 mesh/shape.py"));
#endif
    Shape::allocate();
}

void Shape::create_RF_coil(const double R, const double H, const int N_legs){
    // 
    check_error(this->is_allocated, "shape is already allocated!");
    check_error(!this->is_set, "shape parameters are not set yet!");
    File file;
    file.open("mesh/shape.py", "w");
    file.write("import numpy as np\n");
    file.write("import mesh\n");
    file.write("lambda_ = %21.14E\n", this->lambda);
    file.write("lc = %21.14E\n", this->lc);
    file.write("lines_list = []\n");
    file.write("segments_list = []\n");
    // Begin shape definition
    check_error(N_legs<2, "invalid number of legs!");
    check_error(R<=0.0, "invalid radius!");
    check_error(H<=0.0, "invalid height!");
    file.write("R = %21.14E\n", R);
    file.write("H = %21.14E\n", H);
    file.write("N_legs = %d\n", N_legs);
    file.write("N = int(np.ceil((2.0*np.pi/N_legs)*(R/lambda_)/lc))\n");
    file.write("if N<2:\n");
    file.write("\tN = 2\n");
    file.write("\tpass\n");
    file.write("d_phi = (2.0*np.pi/N_legs)/N\n");
    file.write("for n in range(N_legs):\n");
    file.write("\tphi = n*N*d_phi\n");
    file.write("\tv1 = mesh.Vertex(R*np.cos(phi), R*np.sin(phi), -H/2.0)\n");
    file.write("\tv2 = mesh.Vertex(R*np.cos(phi), R*np.sin(phi), +H/2.0)\n");
    file.write("\tlines_list.append(mesh.Line(v1, v2, True, n+1))\n");
    //
    file.write("\tfor i in range(N):\n");
    file.write("\t\tphi = (n*N+i+0)*d_phi\n");
    file.write("\t\tv1 = mesh.Vertex(R*np.cos(phi), R*np.sin(phi), -H/2.0)\n");
    file.write("\t\tphi = (n*N+i+1)*d_phi\n");
    file.write("\t\tv2 = mesh.Vertex(R*np.cos(phi), R*np.sin(phi), -H/2.0)\n");
    file.write("\t\tsegment = mesh.Segment(v1, v2, False, 0)\n");
    file.write("\t\tsegments_list.append(segment)\n");
    file.write("\t\tcontinue\n");
    //
    file.write("\tfor i in range(N):\n");
    file.write("\t\tphi = (n*N+i+0)*d_phi\n");
    file.write("\t\tv1 = mesh.Vertex(R*np.cos(phi), R*np.sin(phi), +H/2.0)\n");
    file.write("\t\tphi = (n*N+i+1)*d_phi\n");
    file.write("\t\tv2 = mesh.Vertex(R*np.cos(phi), R*np.sin(phi), +H/2.0)\n");
    file.write("\t\tsegment = mesh.Segment(v1, v2, False, 0)\n");
    file.write("\t\tsegments_list.append(segment)\n");
    file.write("\t\tcontinue\n");
    //
    file.write("\tcontinue\n");
    // End of shape definition
    file.write("shape = mesh.Shape(lambda_, lc)\n");
    file.write("for line in lines_list:\n");
    file.write("\tshape.add_line(line)\n");
    file.write("\tcontinue\n");
    file.write("for segment in segments_list:\n");
    file.write("\tshape.add_segment(segment)\n");
    file.write("\tcontinue\n");
    file.write("shape.mesh()\n");
    file.close();
#ifdef _WIN64
    assert(!system("python mesh/shape.py"));
#endif
#ifdef __linux__
    assert(!system("python3 mesh/shape.py"));
#endif
    Shape::allocate();
}