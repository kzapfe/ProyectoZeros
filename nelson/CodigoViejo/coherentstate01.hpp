// Biblioteca que implementa el concepto del estado coherente

#include "simplectic01.hpp"
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <cmath>





# ifndef __coherent__
# define __coherent__

double const pi=3.141592654; //Do you really need an explanation for this?
double const hbarra=0.0750;


class coherent{
  //Clase que implementa Muchas propiedades de un estado
  // coherente.
  // Las convenciones estan en convenciones.pdf

public: 
  
  simplectic zhi;
  //zhi es normalmente el desplazamiento.

  coherent(simplectic &zhi){
    // Un estado coherente centrado en zhi
    this->zhi.q =zhi.q;
    this->zhi.p =zhi.p;
    //cout<<"This is zhi, the Center"<<endl;
    //cout<<zhi.q<<"\t"<<zhi.p<<endl;
  };


  /*
  coherent(simplectic *zhi){
    // Un estado coherente centrado en zhi
    zhi=zhi;
    //cout<<"This is zhi, the Center"<<endl;
    //cout<<zhi.q<<"\t"<<zhi.p<<endl;
  };
  */

  gsl_complex qspace(double x=0.0){
    //La funcion de onda en el espacio q
    // x aqui es una  posicion

    gsl_complex OndaPhi;
    double rho;
    double phase;
    rho=exp(-pow((x-zhi.q),2)/(2.0*hbarra)); 
    phase=(x*zhi.p/(2.0*hbarra));    
    GSL_SET_COMPLEX(&OndaPhi,rho*cos(phase),rho*sin(phase));
    
    return OndaPhi;


  };

  
  gsl_complex pspace(double x=0.00){
    //La funcion de onda en el espacio p
    // x aqui es un momento
    gsl_complex OndaPsi;
    double rho;
    double phase;
    rho=exp(-pow((2.0*x-zhi.p),2)/(4.0*hbarra));
    phase=(2.0*x-zhi.p)*zhi.q/(2.0*hbarra);    
    GSL_SET_COMPLEX(&OndaPsi,rho*cos(phase),rho*sin(phase));
    return OndaPsi;
  };



  double Wigner(simplectic &x){
    //WignerFunction  for the coherent state
    // in 1 degree of freedom  q,p space
    //Definitivamente, nunca entenderas cuando usar los & y los *
    //en las variables. Al parecer es bastante aleetorio.
    
    double result;

    //cout<<x.p-zhi.p<<endl;
    result=1.0/sqrt(pi*hbarra)*
      exp((-1.0/hbarra)*
	       ((x.q-zhi.q)*(x.q-zhi.q)
		+(x.p-zhi.p)*(x.p-zhi.p)));
    return result;
    //Ironically, the simpler one  yet.

  };

  gsl_complex cuerdas(simplectic &x){
    //Aqui x es una Cuerda, whatever that means
    
    double amplitud, fase;
    gsl_complex result;

    amplitud=1.0/(2.0*sqrt(pi*hbarra))
      *exp(-1/(2.0*hbarra)*(x.q*x.q+x.p*x.p));
    //Convencion de Alfredo para el producto simplectico
    fase=1.0/hbarra*(x.q*zhi.p-x.p*zhi.q);    
    GSL_SET_COMPLEX(&result,amplitud*cos(fase),amplitud*sin(fase));
    
    return result;
  };
  
      

};


#endif
