//*Reciclado de las primeras pruebas.

#include <armadillo>
#include "QuantumConstants.hpp"
#include "simplectic01.hpp"


using namespace arma;

//Use your head Karel
//To populate a BALL you dont care if it is symplectic!
//Actually, in your approach,you SCREW simplecticity!!
rowvec CreateRandomPointinBall(int nivel=0, int dimension=2){
  /*Montecarleando la capa de energia numero nivel
  //Mete uno con distribucion uniforme entre dos capas de energia
  Y distribucion Uniforme en el angulo */

  /*Convencion es H=\sum_{dimension } 1/2(x^2+p^2)
   m=1, omega=1 */
  double energia=hbar*((double)nivel+0.5*dimension);
  double energiaaux=hbar*((double)(nivel+1)+0.5*dimension);
  //Marsaglia algoritm
  rowvec x(dimension);
  double radio, angulo;
  double centro, centroaux; 
  double delta;
  
  centro=sqrt(2.*energia);
  centroaux=sqrt(2.*energiaaux);
  delta=centro-centroaux;

  x.randn();
  x=x/norm(x,2);
  // aqui ya tenemos un punto en la esfera de radio 1
  //El radio real es aproximadamente sqrt(2*E)

  radio=as_scalar(randu(1))*delta;
  radio=pow(radio+centro, 1.0/(double)dimension);
  radio=radio*pow(delta+centro, (double)(dimension-1)/(double)dimension);
  
  x=x*radio;
  
  return x;
  
};


mat Populate2BallShell(int nivel=0, int muestreo=1){
  //Crea un array del tamano de muestreo
  // poblando gausianamente cerca de la capa de energia que
  // nos interesa, dada por nivel.
  
  mat points(muestreo, 2);
  
  for(int i=0; i<muestreo; i++) {
    points.row(i)= CreateRandomPointinBall(nivel, 2);
  };
  cout<<"marico"<<endl;
  cout<<points<<endl;
  return points;
    
};


mat Populate4BallShell(int nivel=0, int muestreo=1){
  //Crea un array del tamano de muestreo
  // poblando gausianamente cerca de la capa de energia que
  // nos interesa, dada por nivel.
  
  mat points(muestreo, 4);
  
  for(int i=0; i<muestreo; i++) {
    points.row(i)= CreateRandomPointinBall(nivel, 4);
  };
  //  cout<<"marico"<<endl;
  //cout<<points<<endl;
  return points;
    
};


int signum(double x){
  int sign;
  if(x>0){
    sign=1;
  }else if(x<0){
    sign=-1;
  }else{
    sign=0;
      };
  
  return sign;  
};


