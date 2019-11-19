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

#include "Espacio.h"

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

using namespace std;

int main(int argc, char *argv[])
{
    using clk = chrono::high_resolution_clock;
    auto t1 = clk :: now();

    try {
        if(argc != 5) throw invalid_argument("Wrong arguments \nCorrect use: \nnasteroids-seq num_asteroides num_iteraciones num_planetas semilla" );

        for(int i=0; i<argc ; ++i)
        {
            if(atoi(argv[i]) < 0) throw invalid_argument("Arguments cannot be negative");
        }
    }

    catch(invalid_argument &e)
    {
        cerr << "Exception: " << e.what() << endl;
        return -1;
    }

    num_asteroides = atoi(argv[1]); num_iteraciones = atoi(argv[2]); num_planetas = atoi(argv[3]); seed = atoi(argv[4]);
    
    Espacio main(num_asteroides, num_iteraciones, num_planetas, seed);   

    for(int k = 0; k < num_iteraciones; ++k) main.simulacion();

    auto t2 = clk :: now();
    auto diff = chrono::duration_cast<chrono::microseconds>(t2-t1);
    //cout << "\nTime = " << diff.count()/1000000.0 << " segundos" << endl;
    ofstream coso;
    coso.open("test.txt", ios_base::app);
    coso << diff.count()/1000000.0 << endl;
    coso.close();

    return 0;
}

/**/