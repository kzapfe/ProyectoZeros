#include <armadillo>
#include "simplectic01.hpp"


using namespace arma;

//Constantes del programa
// y del Universo.
double const pi=3.14151926535898; 
double const hbarra=1.000;

simplectic CreateRandomHighEnergyPoint(int nivel=0){
  //Montecarleando la capa de energia numero nivel
  //Mete uno con distribucion Gaussiana en la energia
  //Y distribucion Uniforme en el angulo
  simplectic x;
  double radio, angulo;
  double centro; 
  double delta;
  centro=sqrt(hbarra * nivel);
  delta=sqrt(hbarra * (nivel+1))-sqrt(hbarra * nivel);
 
  radio=as_scalar(randu(1));
  radio=radio*delta+centro;
  
  angulo=2.0*pi*as_scalar(randu(1));

  x.q=cos(angulo)*radio;
  x.p=sin(angulo)*radio;

  return x;


};


simplectic  *PopulateEnergyLevel(int nivel=0, int muestreo=1){
  //Crea un array del tamano de muestreo
  // poblando gausianamente cerca de la capa de energia que
  // nos interesa, dada por nivel.
  
  simplectic *x;
  x=new simplectic[muestreo];

  for(int i=0; i<muestreo; i++) {
      x[i]= CreateRandomHighEnergyPoint(nivel);
  };

  return x;
    
};
