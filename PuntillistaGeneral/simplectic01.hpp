//# include <cmath> Si no lo usas, no lo metas.

//#include <gsl/gsl_complex.h>
//#include <gsl/gsl_complex_math.h>

# ifndef __simplectic__
# define __simplectic__


using namespace std;

#include <algorithm> // este es necesario para el metodo swap

class simplectic{
// a simplectic pair

 
public:
  double q,p; //Notacion usual sentido usual.


   simplectic(double q=0,double p=0) {
     /*Le apunta a los mismas variables, de forma que 
      toda las funciones esten referidas a la 
      variable y NO a su valor */
   this->q =q; this->p =p;  };

  simplectic(const simplectic &x) {
    /*Este segundo constructor le apunta 
      a otro vector, de forma que podamos inicializar
      un vector pasandole los valores de otro*/
    simplectic v;
    q=v.q; p=v.p; };
  

  simplectic(simplectic & Original){
    //Constructor de Copias.
    this->q=Original.q;
    this->p=Original.p;
  };



  friend void swap(simplectic& Uno, simplectic & Dos){
    //intercambia dos CoherentStates, depende de <algorithm>
    using std::swap;
    swap(Uno.q, Dos.q);    
    swap(Uno.p, Dos.p);   
  };

  simplectic& operator=(simplectic RightHandSide){
    swap(*this, RightHandSide);
    return *this;
  };

  simplectic operator+(const simplectic& sumando){
    simplectic result;
    result.q=q+sumando.q;
    result.p=p+sumando.p;
    return result;

  };


  simplectic operator-(const simplectic& restando){
    simplectic result;
    result.q=q-restando.q;
    result.p=p-restando.p;
    return result;

  };


  void SetQandP(double Qaux, double Paux){
    q=Qaux;
    p=Paux;
    return;

  };

  /* avanzar en el vacio es una propiedad
     de la particula, la patada es una propiedad
     del sistema */
  
  void avanzarvacio(double tau=1.00){
    /*avanza tau unidades en un espacio vacio.
      Es una transformacion de deslape. 
      Solo funciona en cartesianas, pendejo.*/ 
    q=q+tau*p;
   

    return;
  };

  
  double simplecticproduct(simplectic & chi){
    //convencion de Alfredo, primero p luego q
    double prod=p*chi.q-q*chi.p;
    return prod;

  };
  



  
  /*  

  gsl_complex ondaplana(double energia=1.0, double hplanck=1.0, double t=0)
  { //Onda plana con valores q,p, E, y hbarra
  
    gsl_complex II;
    GSL_SET_COMPLEX(&II, 0,1);
    gsl_complex auxiliar;
    GSL_SET_COMPLEX(&auxiliar,(simplecticproduct()-energia*t)/hplanck,0);

    gsl_complex ondaplana;
    ondaplana=gsl_complex_exp(gsl_complex_mul(II,auxiliar));
   
       return ondaplana;
    
  };

  */

};

#endif 
