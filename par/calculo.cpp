#include <omp.h>

#include "Espacio.h"

#define G 6.674e-5
#define T 0.1
#define DMIN 5
#define WIDTH 200
#define HEIGHT 200
#define MASA 1000
#define SDM 50

void Espacio::simulacion()
{
    vector<double> fuerza_x(num_asteroides, 0.0);
	vector<double> fuerza_y(num_asteroides, 0.0); 

    calcular_fuerzas(fuerza_x, fuerza_y);
    calcular_posicion(fuerza_x, fuerza_y);
    calculo_rebotes();
    print_end();
}

Espacio::Espacio(int na, int ni, int np, int s) :
    asteroides{}, planetas{},
    num_asteroides{na}, num_iteraciones{ni}, num_planetas{np}, seed{s}
{
    ofstream init_conf("init_conf.txt");
    init_conf.setf(ios_base::fixed, ios_base::floatfield);
    init_conf.precision(3);

    if (!init_conf.is_open())
    {
        cout << "Unable to open init_conf.txt" << "\n";
        exit(-1);
    }

    init_conf << num_asteroides << " " << num_iteraciones << " " << num_planetas << " " << seed << endl;
    
    asteroides.reserve(num_asteroides);
	planetas.reserve(num_planetas); 
    
    // Random distributions
    default_random_engine re{seed};
    uniform_real_distribution<double> xdist{0.0, nextafter(WIDTH, numeric_limits<double>::max())};
    uniform_real_distribution<double> ydist{0.0, nextafter(HEIGHT, numeric_limits<double>::max())};
    normal_distribution<double> mdist{MASA, SDM};

    for(int i = 0 ; i < num_asteroides ; ++i)
    {
        Asteroide aux = Asteroide{xdist(re), ydist(re), mdist(re)};
        init_conf << aux.posx << " " << aux.posy << " " << aux.masa << endl;

        asteroides.push_back(aux);
    }

    for(int i = 0 ; i < num_planetas ; ++i)
    {
        Planeta aux;

        switch(i % 4)
        {
            case 0: aux = Planeta{0.0, ydist(re), mdist(re)*10}; break;
            case 1: aux = Planeta{xdist(re), 0.0, mdist(re)*10}; break;
            case 2: aux = Planeta{HEIGHT, ydist(re), mdist(re)*10}; break;
            case 3: aux = Planeta{xdist(re), WIDTH, mdist(re)*10}; break;
        }
        init_conf << aux.posx << " " << aux.posy << " " << aux.masa << endl;

        planetas.push_back(aux);
    }
    init_conf.close();
}

double Espacio::angulo(double xdif, double ydif)
{
    double ang; 
    // Cálculo de la pendiente
    double m = ydif / xdif;

    if(m < -1) m = -1; // Si es menor que -1 
    if(m > 1) m = 1; // Si es mayor que 1

    // Cálculo del ángulo de influencia
    return ang = atan(m);
}

void Espacio::calcular_fuerzas(vector<double> &fuerza_x, vector<double> &fuerza_y)
{
    #pragma omp parallel for
    for(int i = 0 ; i < num_asteroides ; ++i)
    {
        // Atracción con otros asteroides
        for(int j = i+1 ; j < num_asteroides ; ++j)
        {
            //cout << "\nA " << i << " y A " << j << "\n";
            double dist, ang, fx, fy, xdif, ydif;
            // Cálculo de la distancia
            xdif = asteroides[i].posx - asteroides[j].posx;
            ydif = asteroides[i].posy - asteroides[j].posy; 
            dist = sqrt(pow(xdif,2) + pow(ydif,2));
            if(dist > DMIN)
            {
                // Cálculo del ángulo de influencia
                ang = angulo(xdif, ydif);
                //cout << "Ángulo " << ang << "\n";
                // Cálculo de la fuerza
                double op1 = G * asteroides[i].masa * asteroides[j].masa;
                double op2 = pow(dist, 2);
                fx = (op1/op2)*cos(ang);
                fy = (op1/op2)*sin(ang);
                //cout << "Fuerza " << sqrt(pow(fx,2) + pow(fy,2)) << "\n";
                // Añade las fuerzas ejercidas del asteroide i a j
                fuerza_x[i] += fx; 
                fuerza_y[i] += fy;
                // Añade las fuerzas ejercidas del asteroide i a j en dirección opuesta
                fuerza_x[j] -= fx; 
                fuerza_y[j] -= fy;
            }
        }
        // Atracción con planetas
        for(int j = 0 ; j < num_planetas ; ++j)
        {
            //cout << "\nA " << i << " y P " << j << "\n";
            double dist, ang, fx, fy, xdif, ydif;
            // Cálculo de la distancia
            xdif = asteroides[i].posx - planetas[j].posx;
            ydif = asteroides[i].posy - planetas[j].posy; 
            dist = sqrt(pow(xdif,2) + pow(ydif,2));
            // Cálculo del ángulo de influencia
            ang = angulo(xdif, ydif);
            //cout << "Ángulo " << ang << "\n";
            // Cálculo de la fuerza
            double op1 = G * asteroides[i].masa * planetas[j].masa;
            double op2 = pow(dist, 2);
            fx = (op1/op2)*cos(ang);
            fy = (op1/op2)*sin(ang);
            //cout << "Fuerza " << sqrt(pow(fx,2) + pow(fy,2)) << "\n";
            // Añade las fuerzas ejercida por el planeta al asteroide
            fuerza_x[i] += fx; 
            fuerza_y[i] += fy;
        }
        if(fuerza_x[i] > 100) fuerza_x[i] = 100;
        if(fuerza_y[i] > 100) fuerza_y[i] = 100;
    }
}

void Espacio::calcular_posicion(vector<double> &fuerza_x, vector<double> &fuerza_y)
{
    #pragma omp parallel for
    for(int i = 0 ; i < num_asteroides ; ++i)
    {
        double ax, ay;

        // Cálculo aceleración
        ax = fuerza_x[i] / asteroides[i].masa; 
        ay = fuerza_y[i] / asteroides[i].masa;
        
        // Cálculo velocidad
        asteroides[i].velx += ceil(ax*T*100000)/100000;
        asteroides[i].vely += ceil(ay*T*100000)/100000;

        // Cálculo posición
        asteroides[i].posx += ceil(asteroides[i].velx*T*100000)/100000;
        asteroides[i].posy += ceil(asteroides[i].vely*T*100000)/100000;

        // Rebote con los bordes
        if(asteroides[i].posx <= 0.0) 
        { 
            asteroides[i].posx = 5.0; 
            asteroides[i].velx = -asteroides[i].velx;
        }
        else if(asteroides[i].posy <= 0.0) 
        { 
            asteroides[i].posy = 5.0; 
            asteroides[i].vely = -asteroides[i].vely;  
        }
        else if(asteroides[i].posx >= WIDTH) 
        { 
            asteroides[i].posx = WIDTH - 5.0; 
            asteroides[i].velx = -asteroides[i].velx;
        }
        else if(asteroides[i].posy >= HEIGHT) 
        { 
            asteroides[i].posy = HEIGHT - 5.0; 
            asteroides[i].vely = -asteroides[i].vely; 
        }
    }
}

void Espacio::calculo_rebotes()
{
    for(int i = 0 ; i < num_asteroides ; ++i)
    {
        for(int j = i+1 ; j < num_asteroides ; ++j)
        {
            double dist, xdif, ydif;

            // Cálculo de la distancia
            xdif = asteroides[i].posx - asteroides[j].posx;
            ydif = asteroides[i].posy - asteroides[j].posy; 
            dist = sqrt(pow(xdif,2) + pow(ydif,2));

            if(dist < DMIN)
            {
                // Variables auxiliares
                double vx = asteroides[j].velx; 
                double vy = asteroides[j].vely;
                // Se intercambian las velocidades j <- i
                asteroides[j].velx = asteroides[i].velx; 
                asteroides[j].vely = asteroides[i].vely;
                // Se intercambian las velocidades i <- j
                asteroides[i].velx = vx; 
                asteroides[i].vely = vy;
            } 
        }
    }
}

void Espacio::print_end()
{
    std::ofstream out("out.txt"); 
    out.setf(std::ios_base::fixed, std::ios_base::floatfield);
    out.precision(3);

    for(int i = 0 ; i < num_asteroides ; ++i)
        out << asteroides[i].posx << " " << asteroides[i].posy << " " << asteroides[i].velx << " " 
           << asteroides[i].vely << " " << asteroides[i].masa << " " << std::endl;
    out.close();
}