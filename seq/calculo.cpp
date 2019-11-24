#include "Espacio.h"

#define G 6.674e-5
#define T 0.1
#define DMIN 5
#define WIDTH 200
#define HEIGHT 200
#define MASA 1000
#define SDM 50

using namespace std;

/**
 * Constructor de la clase Espacio:
 *  - Guarda los parámetros de entrada
 *  - Inicializa los vectores de asteroides y planetas
 *  - Escribe el archivo inicial
 *  - Crea los asteroides y planetas, y los guarda en su vector correspodiente
 * 
 * @param na Número de asteroides
 * @param ni Número de iteraciones
 * @param np Número de planetas
 * @param s  Semilla de la distribución aleatoria
 */
Espacio::Espacio(int na, int ni, int np, int s) :
    asteroides{}, planetas{},
    num_asteroides{na}, num_iteraciones{ni}, num_planetas{np}, seed{s}
{
    // Abre el archivo inicial con precisión 3
    ofstream init_conf("init_conf.txt");
    init_conf.setf(ios_base::fixed, ios_base::floatfield);
    init_conf.precision(3);
    // Comprueba que el archivo está abierto
    if (!init_conf.is_open())
    {
        // Devuelve error en caso contrario
        cout << "Unable to open init_conf.txt" << "\n";
        exit(-1);
    }
    // Escribe los parámetros de entrada en el archivo inicial
    init_conf << num_asteroides << " " << num_iteraciones << " " << num_planetas << " " << seed << endl;
    // Reserva el espacio de los vectores de asteroides y planetas
    asteroides.reserve(num_asteroides);
	planetas.reserve(num_planetas); 
    
    // Distribución aleatoria
    default_random_engine re{seed};
    uniform_real_distribution<double> xdist{0.0, nextafter(WIDTH, numeric_limits<double>::max())};
    uniform_real_distribution<double> ydist{0.0, nextafter(HEIGHT, numeric_limits<double>::max())};
    normal_distribution<double> mdist{MASA, SDM};
    // Crea los asteroides con sus parámetros correspondientes y los guarda en su vector
    for(int i = 0 ; i < num_asteroides ; ++i)
    { // Start loop
        // Llamada al constructor de la clase Asteroide. Recibe como parámetros las posiciones x e y yy la masa
        Asteroide ast = Asteroide{xdist(re), ydist(re), mdist(re)};
        // Escribe los datos del asteroide en el archivo inicial
        init_conf << ast.posx << " " << ast.posy << " " << ast.masa << endl;
        // Añade el nuevo asteroide al final del vector
        asteroides.push_back(ast);
    } // End loop

    // Crea los planetas con sus parámetros correspodientes y los guarda en su vector
    for(int i = 0 ; i < num_planetas ; ++i)
    {// Start loop
        Planeta pla;
        // Switch para distribuir los planetas por los bordes del espacio
        switch(i % 4)
        {   
            // Planeta en el borde izquierdo
            case 0: aux = Planeta{0.0, ydist(re), mdist(re)*10}; break;
            // Planeta en el borde superior
            case 1: aux = Planeta{xdist(re), 0.0, mdist(re)*10}; break;
            // Planeta en el borde izquierdo
            case 2: aux = Planeta{HEIGHT, ydist(re), mdist(re)*10}; break;
            // Planeta en el borde inferior
            case 3: aux = Planeta{xdist(re), WIDTH, mdist(re)*10}; break;
        }
        // Escribe los planetas en el archivo inicial
        init_conf << pla.posx << " " << pla.posy << " " << pla.masa << endl;
        // Añade el nuevo planeta al final de su vector
        planetas.push_back(pla);
    } // End loop
    // Cierra el fichero
    init_conf.close();
}

/**
 * Inicia la simulación del movimiento de los asteroides
 */
void Espacio::simulacion()
{
    // Crea los vectores de las fuerzas en cada eje y los inicializa en 0
    vector<double> fuerza_x(num_asteroides, 0.0);
	vector<double> fuerza_y(num_asteroides, 0.0); 
    // Llamada a las funciones de la clase Espacio
    calcular_fuerzas(fuerza_x, fuerza_y);
    calcular_posicion(fuerza_x, fuerza_y);
    calculo_rebotes();
    print_end();
}

/**
 * Calcula la fuerza ejercida sobre los asteroides por otros asteroides y por los planetas. Estas fuerzas
 * son guardadas en los vectores en su eje correspondiente
 * @param fuerza_x Vector con las fuerzas en el eje x de los asteroides
 * @param fuerza_y Vector con las fuerzas en el eje y de los asteroides
 */
void Espacio::calcular_fuerzas(vector<double> &fuerza_x, vector<double> &fuerza_y)
{
    for(int i = 0 ; i < num_asteroides ; ++i)
    { // Start loop 1
        // Atracción con otros asteroides 
        for(int j = i+1 ; j < num_asteroides ; ++j)
        { // Start loop 2
            double dist, ang, fx, fy, xdif, ydif;
            // Cálculo de la distancia
            xdif = asteroides[i].posx - asteroides[j].posx;
            ydif = asteroides[i].posy - asteroides[j].posy; 
            dist = sqrt(pow(xdif,2) + pow(ydif,2));
            // Comprueba que la distancia entre los asteroides es la permitida para no rebotar
            if(dist > DMIN)
            {
                // Cálculo del ángulo de influencia
                ang = angulo(xdif, ydif);
                // Cálculo de la fuerza
                double op1 = G * asteroides[i].masa * asteroides[j].masa;
                double op2 = pow(dist, 2);
                fx = (op1/op2)*cos(ang);
                fy = (op1/op2)*sin(ang);
                // Añade a los vectores las fuerzas ejercidas del asteroide i a j
                fuerza_x[i] += fx; 
                fuerza_y[i] += fy;
                // Añade a los vectores las fuerzas ejercidas del asteroide i a j en dirección opuesta
                fuerza_x[j] -= fx; 
                fuerza_y[j] -= fy;
            }
        } // End loop 2
        // Atracción con planetas
        for(int j = 0 ; j < num_planetas ; ++j)
        { // Start loop 3
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
            // Añade a los vectores las fuerzas ejercidas por el planeta al asteroide
            fuerza_x[i] += fx; 
            fuerza_y[i] += fy;
        } // End loop 3
        // Comprueba que la fuerza total ejercida sobre un asteroide no supere la fuerza máxima permitida
        if(fuerza_x[i] > 100) fuerza_x[i] = 100;
        if(fuerza_y[i] > 100) fuerza_y[i] = 100;
    } // Start loop 1
}

/**
 * Se dedica a calcular la posición de los asteroides teniendo en cuenta su velocidad
 * @param fuerza_x Vector con las fuerzas en el eje x de los asteroides
 * @param fuerza_y Vector con las fuerzas en el eje y de los asteroides
 */
void Espacio::calcular_posicion(vector<double> &fuerza_x, vector<double> &fuerza_y)
{
    for(int i = 0 ; i < num_asteroides ; ++i)
    { // Start loop
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
        // Rebote con los bordes: reajusta la posición de los asteroides en caso de que estén fuera de los márgenes
        // especificados por el problema
        if(asteroides[i].posx <= 0.0) // Si traspasa el borde izquierdo
        { 
            asteroides[i].posx = 5.0; 
            asteroides[i].velx = -asteroides[i].velx;
        }
        else if(asteroides[i].posy <= 0.0) // Si traspasa el borde superior
        { 
            asteroides[i].posy = 5.0; 
            asteroides[i].vely = -asteroides[i].vely;  
        }
        else if(asteroides[i].posx >= WIDTH) // Si traspasa el borde derecho
        { 
            asteroides[i].posx = WIDTH - 5.0; 
            asteroides[i].velx = -asteroides[i].velx;
        }
        else if(asteroides[i].posy >= HEIGHT) // Si traspasa el borde inferior
        { 
            asteroides[i].posy = HEIGHT - 5.0; 
            asteroides[i].vely = -asteroides[i].vely; 
        }
    } // End loop
}

/**
 * Evalúa el posible rebote de los asteroides. Este rebote ocurre cuando dos asteroides se encuentran a menos de 
 * distancia 5. En caso de que así sea, estos intercambian sus velocidades.
 */
void Espacio::calculo_rebotes()
{
    for(int i = 0 ; i < num_asteroides ; ++i)
    { // Start loop 1
        for(int j = i+1 ; j < num_asteroides ; ++j)
        { // Start loop 2
            double dist, xdif, ydif;
            // Cálculo de la distancia
            xdif = asteroides[i].posx - asteroides[j].posx;
            ydif = asteroides[i].posy - asteroides[j].posy; 
            dist = sqrt(pow(xdif,2) + pow(ydif,2));
            // Comprueba que la distancia entre los asteroides es la mínima para chocar
            if(dist <= DMIN)
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
        } // End loop 2
    } // End loop 1
}

/**
 * Calcula el ángulo de influencia de la fuerza de un asteroide o planeta sobre un asteroide
 * @param xdif Diferencia en la distancia en el eje x
 * @param ydif Diferencia en la distancia en el eje y
 */
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
 * Imprime el archivo de salida del programa el cual recoge las posiciones, velocidades y masa finales de los asteroides
 */
void Espacio::print_end()
{
    // Crea el archivo de salida del programa
    ofstream out("out.txt"); 
    out.setf(ios_base::fixed, ios_base::floatfield);
    out.precision(3);
    // Comprueba que el archivo esté abierto, en caso incorrecto error y fin de la ejecución
    if (!out.is_open())
    {
        cout << "Unable to open out.txt" << "\n";
        exit(-1);
    }
    // Ecribe las posiciones, velocidades y masa de los asteroides en el archivo de salida
    for(int i = 0 ; i < num_asteroides ; ++i)
        out << asteroides[i].posx << " " << asteroides[i].posy << " " << asteroides[i].velx << " " 
           << asteroides[i].vely << " " << asteroides[i].masa << " " << endl;
    // Cierra el archivo
    out.close();
}