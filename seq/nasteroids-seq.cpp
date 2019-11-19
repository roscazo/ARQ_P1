/**
 * @authors Rodrigo Morales Sánchez & Marina García Caro
 * Arquitectura de Computadoras - Práctica 1
 * Oct/Nov 2019
 **/

#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <stdexcept>
#include <fstream>
#include <chrono>

#include "Elemento.h"

#define G 6.674e-5
#define T 0.1
#define DMIN 5
#define WIDTH 200
#define HEIGHT 200
#define MASA 1000
#define SDM 50

int num_asteroides;
int num_iteraciones;
int num_planetas;
int seed;

double angulo(double xdif, double ydif)
{
    double ang; 
    // Cálculo de la pendiente
    double m = ydif / xdif;

    if(m < -1) m = -1; // Si es menor que -1 
    if(m > 1) m = 1; // Si es mayor que 1

    // Cálculo del ángulo de influencia
    return ang = atan(m);
}

void calcular_fuerzas(vector<double> &fuerza_x, vector<double> &fuerza_y)
{
    for(int i = 0 ; i < num_asteroides ; ++i)
    {
        // Atracción con otros asteroides
        for(int j = i+1 ; j < num_asteroides ; ++j)
        {
            double dist, ang, fx, fy, xdif, ydif;
            // Cálculo de la distancia
            xdif = asteroides[i].posx - asteroides[j].posx;
            ydif = asteroides[i].posy - asteroides[j].posy; 
            dist = sqrt(pow(xdif,2) + pow(ydif,2));

            if(dist > DMIN)
            {
                // Cálculo del ángulo de influencia
                ang = angulo(xdif, ydif);

                // Cálculo de la fuerza
                double op1 = G * asteroides[i].masa * asteroides[j].masa;
                double op2 = pow(dist, 2);
                fx = (op1/op2)*cos(ang);
                fy = (op1/op2)*sin(ang);

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
            double dist, ang, fx, fy, xdif, ydif;
            // Cálculo de la distancia
            xdif = asteroides[i].posx - asteroides[j].posx;
            ydif = asteroides[i].posy - asteroides[j].posy; 
            dist = sqrt(pow(xdif,2) + pow(ydif,2));

            // Cálculo del ángulo de influencia
            ang = angulo(xdif, ydif);
            
            // Cálculo de la fuerza
            double op1 = G * asteroides[i].masa * planetas[j].masa;
            double op2 = pow(dist, 2);
            fx = (op1/op2)*cos(ang);
            fy = (op1/op2)*sin(ang);

            // Añade las fuerzas ejercida por el planeta al asteroide
            fuerza_x[i] += fx; 
            fuerza_y[i] += fy;
        }
    }
    if(fuerza_x[i] > 100) fuerza_x[i] = 100;
    if(fuerza_y[i] > 100) fuerza_y[i] = 100;
}

void calcular_posicion(vector<double> &fuerza_x, vector<double> &fuerza_y)
{
    for(int i = 0 ; i < num_asteroides ; ++i)
    {
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

void calculo_rebotes()
{
    for(int i = 0 ; i < num_asteroides ; ++i)
    {
        for(int j = i+1 ; j < num_asteroides ; ++j)
        {
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

int main(int argc, char *argv[])
{
    using clk = std::chrono::high_resolution_clock;
    auto t1 = clk :: now();

    try {
        if(argc != 5) throw std::invalid_argument("Wrong arguments \nCorrect use: \nnasteroids-seq num_asteroides num_iteraciones num_planetas semilla" );

        for(int i=0; i<argc ; ++i)
        {
            if(atoi(argv[i]) < 0) throw std::invalid_argument("Arguments cannot be negative");
        }
    }

    catch(std::invalid_argument &e)
    {
        std::cerr << "Exception: " << e.what() << std::endl;
        return -1;
    }

    std::ofstream init_conf("init_conf.txt");
    init_conf.setf(std::ios_base::fixed, std::ios_base::floatfield);
    init_conf.precision(3);

    if (!init_conf.is_open())
    {
        std::cout << "Unable to open init_conf.txt" << "\n";
        return -1;
    }

    num_asteroides = atoi(argv[1]); num_iteraciones = atoi(argv[2]); num_planetas = atoi(argv[3]); seed = atoi(argv[4]);
    init_conf << num_asteroides << " " << num_iteraciones << " " << num_planetas << " " << seed << std::endl;
    
    std::vector<Asteroide> asteroides(num_asteroides);
	std::vector<Planeta> planetas(num_planetas); 
    
    // Random distributions
    std::default_random_engine re{seed};
    std::uniform_real_distribution<double> xdist{0.0, nextafter(WIDTH, std::numeric_limits<double>::max())};
    std::uniform_real_distribution<double> ydist{0.0, nextafter(HEIGHT, std::numeric_limits<double>::max())};
    std::normal_distribution<double> mdist{MASA, SDM};

    for(int i = 0 ; i < num_asteroides ; ++i)
    {
        asteroides[i] = Asteroide{xdist(re), ydist(re), mdist(re)};
        init_conf << asteroides[i].posx << " " << asteroides[i].posy << " " << asteroides[i].masa << std::endl;
    }

    for(int i = 0 ; i < num_planetas ; ++i)
    {
        switch(i % 4)
        {
            case 0: planetas[i] = Planeta{0.0, ydist(re), mdist(re)}; break;
            case 1: planetas[i] = Planeta{xdist(re), 0.0, mdist(re)}; break;
            case 2: planetas[i] = Planeta{HEIGHT, ydist(re), mdist(re)}; break;
            case 3: planetas[i] = Planeta{xdist(re), WIDTH, mdist(re)}; break;
        }
        init_conf << planetas[i].posx << " " << planetas[i].posy << " " << planetas[i].masa << std::endl;
    }
    init_conf.close();

    std::vector<double> fuerza_x(num_asteroides, 0.0);
	std::vector<double> fuerza_y(num_asteroides, 0.0);    

    for(int k = 0; k < num_iteraciones; ++k)
    {
        calcular_fuerzas(fuerza_x, fuerza_y);
        calcular_posicion(fuerza_x, fuerza_y);
        calcular_rebotes();
    }

    std::ofstream out("out2.txt"); 
    out.setf(std::ios_base::fixed, std::ios_base::floatfield);
    out.precision(3);

    for(int i = 0 ; i < num_asteroides ; ++i)
        out << asteroides[i].pos_x << " " << asteroides[i].pos_y << " " << asteroides[i].vel_x << " " 
           << asteroides[i].vel_y << " " << asteroides[i].masa << " " << std::endl;
    out.close();

    auto t2 = clk :: now();
    auto diff = std::chrono::duration_cast<std::chrono::microseconds>(t2-t1);
    //std::cout << "\nTime = " << diff.count()/1000000.0 << " segundos" << std::endl;
    std::ofstream coso;
    coso.open("test.txt", std::ios_base::app);
    coso << diff.count()/1000000.0 << std::endl;
    coso.close();

    return 0;
}

/**/