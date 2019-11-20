#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <stdexcept>
#include <fstream>
#include <chrono>

using namespace std;

class Asteroide
{
public:
    double posx;
    double posy; 
    double masa;
    double velx = 0.0;
    double vely = 0.0;

    Asteroide(double x, double y, double m) :
        posx{x}, posy{y}, masa{m} {};

    Asteroide() {};    

    ~Asteroide() {};   

};

class Planeta
{
public:
    double posx;
    double posy; 
    double masa;

    Planeta(double x, double y, double m) :
        posx{x}, posy{y}, masa{m} {};

    Planeta() {};

    ~Planeta() {};

};

class Espacio
{
public:

    void simulacion();
    Espacio(int na, int ni, int np, int s);
    ~Espacio() {};

private:
    vector<Asteroide> asteroides;
	vector<Planeta> planetas;

    int num_asteroides;
    int num_iteraciones;
    int num_planetas;
    int seed;   

    double angulo(double xdif, double ydif);
    void calcular_fuerzas(vector<double> &fuerza_x, vector<double> &fuerza_y);
    void calcular_posicion(vector<double> &fuerza_x, vector<double> &fuerza_y);
    void calculo_rebotes();
    void print_end();

};