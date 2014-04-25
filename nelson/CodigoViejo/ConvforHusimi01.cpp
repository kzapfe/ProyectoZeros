/*
Este programa saca la 
Distribucion de Hussimi a partir de una distribucion de Wigner.

Esto es, obtiene la convolucion

 */
// estos dos son los de gsl para rutinas aleatoreas
//#include <gsl_rng.h>
//#include <gsl_randist.h>
// De momento usaremos solo la guassiana de armadillo

#include <cstdlib>
#include <armadillo>
#include "simplectic01.hpp"
#include "ParametrosGlobales.hpp"
#include "RutinasNelson03.hpp"


using namespace std;



int main(){


  

  //La funcion de Wigner, completa con sus valores
  mat WignerFunction;
  
  WignerFunction.load("sinespacios.dat", raw_ascii );

  // la longitud de los datos originales.
  int longdat=WignerFunction.n_rows;

  double deltaxq, deltaxp;

  deltaxq=WignerFunction(1,0)-WignerFunction(0,0);
  deltaxp=WignerFunction(1,1)-WignerFunction(0,1);
  
  double Husimi=0.00;
  double Xhusq=0.00, Xhusp=0.00;

  cout<<"Tenemos este chichiputal de datos "<< longdat<<endl;
  cout<<"Delta xq="<< deltaxq<<endl;
  cout<<"Delta xp="<< deltaxp<<endl;



  for(int i=0; i<longdat; i++){
    Husimi+=WignerFunction(i,2)*
      exp((-(WignerFunction(i,0)-Xhusq)*(WignerFunction(i,0)-Xhusq)
	   -(WignerFunction(i,1)-Xhusp)*(WignerFunction(i,1)-Xhusp)/hbarra));
  
  };
  
  Husimi=Husimi/(hbarra*pi)*deltaxp*deltaxp;
  
  cout<<"el valor de la funcion de Husimi en el origen es "<<
    Husimi<<endl;



}

