/* Objeto intervalito 
 este objeto clasifica un intervalito de acuerdo a tres
propiedades. dos reales (sus limites) y si una funcion cambia
de signo en estos dos limites
*/



#ifndef __intervalito__
#define  __intervalito__




using namespace std;



class intervalito{
// Un intervalito con un posible cambio de signo
  // de una funcion R->C en el.

public:
  double menor,mayor; //Notacion usual sentido usual.
  bool cambia;

  intervalito(double mayor=0.0,double menor=1.0) {
    /*El constructor por omision produce el intervalo unidad
      variable y NO a su valor La funcion por omision es la identidad
    y el signo de cero es cero.*/
   this->mayor =mayor; this->menor =menor;  
   this->cambia=true;   
 
  };

  intervalito(const intervalito &INTERVAL) {
    
    mayor=INTERVAL.mayor; 
    menor=INTERVAL.menor;  
    cambia=INTERVAL.cambia;   
 
    };

  //  ~intervalito();

};
  
#endif


intervalito *Particion(int partes, intervalito original){
    //Parte el intervalito en  partes intervalitos iguales
    intervalito *pieza;
    double longitud;
    longitud=original.mayor-original.menor;


    pieza=new intervalito[partes];
    
    
    for(int i=0; i<partes; i++){
      double minima, maxima;
      minima=(double)i/double(partes)*(longitud)+original.menor;;
      maxima=(double)(i+1)/double(partes)*(longitud)+original.menor;
      pieza[i].menor=minima;
      pieza[i].mayor=maxima;    
    };
    
    return pieza;

};
