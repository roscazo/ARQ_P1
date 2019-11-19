class Elemento
{
public:
    double posx;
    double posy; 
    double masa;

    Elemento(double x, double y, double m):
        posx{x}, posy{y}, masa{m} {}

    Elemento() = default;

    ~Elemento();

    // virtual void setPos(double x, double y){ posx = x; posy = y; }
    // virtual void setMasa(double m){ masa = m; }

    // double getPos_X(){ return posx; }
    // double getPos_Y(){ return posy; }
    // double getMasa(){ return masa; }
};

class Asteroide : public Elemento
{
public:
    double velx = 0.0;
    double vely = 0.0;

    Asteroide(double x, double y, double m) : Elemento(x, y, m) {};

    Asteroide() = default;   

    // void setVel(double x, double y){ vx = x; vy = y; }
    
    // double getVX(){ return vx; }
    // double getVY(){ return vy; }
};

class Planeta : public Elemento
{
public:
    Planeta(double x, double y, double m) : Elemento(x, y, m*10){};

    Planeta() = default;
};


Elemento::~Elemento()
{
}
