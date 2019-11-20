#include <omp.h>

#include "Espacio.h"

#define G 6.674e-5
#define T 0.1
#define DMIN 5
#define WIDTH 200
#define HEIGHT 200
#define MASA 1000
#define SDM 50

/**
 * Inicia la simulación 
 **/
void Espacio::simulacion()
{
    // Inicializa en 0 los vectores de fuerzas de atracción de los asteroides
    vector<double> fuerza_x(num_asteroides, 0.0);
	vector<double> fuerza_y(num_asteroides, 0.0); 
    // Llamada a las funciones de cálculo
    calcular_fuerzas(fuerza_x, fuerza_y);
    calcular_posicion(fuerza_x, fuerza_y);
    calculo_rebotes();
    print_end();
}

/**
 * Constructor de la clase Espacio
 * Inicializa las variables de entrada y los vectores de asteroides y planetas
 * Además, llena los atributos de estos elementos y los introduce en su vector correspondiente
 **/
Espacio::Espacio(int na, int ni, int np, int s) :
    asteroides{}, planetas{},
    num_asteroides{na}, num_iteraciones{ni}, num_planetas{np}, seed{s}
{
    ofstream init_conf("init_conf.txt");
    init_conf.setf(ios_base::fixed, ios_base::floatfield);
    init_conf.precision(3);
    // Comprueba si el archivo inicial está abierto. En caso contrario, error
    if (!init_conf.is_open())
    {
        cout << "Unable to open init_conf.txt" << "\n";
        exit(-1);
    }
    // Escribe los atributos de entrada en el archivo inicial
    init_conf << num_asteroides << " " << num_iteraciones << " " << num_planetas << " " << seed << endl;
    // Reserva el espacio de los vectores de asteroides y planetas
    asteroides.reserve(num_asteroides);
	planetas.reserve(num_planetas);
    
    // Distribución random
    default_random_engine re{seed};
    uniform_real_distribution<double> xdist{0.0, nextafter(WIDTH, numeric_limits<double>::max())};
    uniform_real_distribution<double> ydist{0.0, nextafter(HEIGHT, numeric_limits<double>::max())};
    normal_distribution<double> mdist{MASA, SDM};
    // Llena el vector asteroides con el número de asteroides de entrada
    for(int i = 0 ; i < num_asteroides ; ++i)
    {
        // Llamada al constructor de Asteroide
        Asteroide aux = Asteroide{xdist(re), ydist(re), mdist(re)};
        // Escribe los atributos del asteroide en el archivo inicial
        init_conf << aux.posx << " " << aux.posy << " " << aux.masa << endl;
        // Push para llenar el vector asteroides
        asteroides.push_back(aux);
    }
    // Llena el vector planetas con el número de planetas de entrada
    for(int i = 0 ; i < num_planetas ; ++i)
    {
        Planeta aux;

        switch(i % 4)
        {
            // Llamada al constructor de Planeta
            case 0: aux = Planeta{0.0, ydist(re), mdist(re)*10}; break;
            case 1: aux = Planeta{xdist(re), 0.0, mdist(re)*10}; break;
            case 2: aux = Planeta{HEIGHT, ydist(re), mdist(re)*10}; break;
            case 3: aux = Planeta{xdist(re), WIDTH, mdist(re)*10}; break;
        }
        // Escribe los atributos del planeta en el archivo inicial
        init_conf << aux.posx << " " << aux.posy << " " << aux.masa << endl;
        // Push para llenar el vector planetas
        planetas.push_back(aux);
    }
    // Cierra el archivo inicial
    init_conf.close();
}

/**
 * Función auxiliar para calcular el ángulo de influencia de una fuerza
 **/
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

/**
 * 
 * 
 **/
void Espacio::calcular_fuerzas(vector<double> &fuerza_x, vector<double> &fuerza_y)
{
    vector<double> x_privada;
	vector<double> y_privada;

    omp_set_num_threads(16);	
	#pragma omp parallel 
    {
        //Var para identificar cada thread
		const int num_threads = omp_get_num_threads();
		const int thread_id = omp_get_thread_num();

		// Creación de vectores privados en cada thread para evitar conflictos y reducción del overhead
		#pragma omp single
		{
			x_privada = vector<double> (num_asteroides*num_threads, 0.0);
			y_privada = vector<double> (num_asteroides*num_threads, 0.0);
		}

        #pragma omp for 
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
                    x_privada[thread_id*num_asteroides+i] += fx; 
                    y_privada[thread_id*num_asteroides+i] += fy;
                    // Añade las fuerzas ejercidas del asteroide i a j en dirección opuesta
                    x_privada[thread_id*num_asteroides+j] -= fx; 
                    y_privada[thread_id*num_asteroides+j] -= fy;
                }
            }

            // Atracción con planetas
            for(int j = 0 ; j < num_planetas ; ++j)
            {
                double dist, ang, fx, fy, xdif, ydif;
                // Cálculo de la distancia
                xdif = asteroides[i].posx - planetas[j].posx;
                ydif = asteroides[i].posy - planetas[j].posy; 
                dist = sqrt(pow(xdif,2) + pow(ydif,2));
                // Cálculo del ángulo de influencia
                ang = angulo(xdif, ydif);
                // Cálculo de la fuerza
                double op1 = G * asteroides[i].masa * planetas[j].masa;
                double op2 = pow(dist, 2);
                fx = (op1/op2)*cos(ang);
                fy = (op1/op2)*sin(ang);
                // Añade las fuerzas ejercida por el planeta al asteroide
                x_privada[thread_id*num_asteroides+i] += fx; 
                y_privada[thread_id*num_asteroides+i] += fy;
            }
        }

        //Reducción manual de los vectores privados
		#pragma omp for 
		for(int l = 0; l < num_asteroides; ++l)
        {
			for(int t = 0; t < num_threads; ++t)
            {
				fuerza_x[l] += x_privada[num_asteroides*t+l];
				fuerza_y[l] += y_privada[num_asteroides*t+l];
                if(fuerza_x[l] > 100) fuerza_x[l] = 100;
                if(fuerza_y[l] > 100) fuerza_y[l] = 100;
			}
		}
    }
}

/**
 * 
 * 
 **/
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
        asteroides[i].velx += ax*T;
        asteroides[i].vely += ay*T;

        // Cálculo posición
        asteroides[i].posx += asteroides[i].velx*T;
        asteroides[i].posy += asteroides[i].vely*T;

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

/**
 * 
 * 
 **/
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

/**
 * 
 * 
 **/
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