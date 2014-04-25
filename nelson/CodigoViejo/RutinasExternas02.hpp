#include <armadillo>
#include "simplectic01.hpp"


using namespace arma;

//Constantes del programa
// y del Universo.
double const pi=3.14151926535898; 
double const hbarra=1.000;


//Use your head Karel
//To populate a BALL you dont care if it is symplectic!
//Actually, in your approach,you SCREW simplecticity!!
rowvec CreateRandomPointinBall(int nivel=0, int dimension=2){
  /*Montecarleando la capa de energia numero nivel
  //Mete uno con distribucion Gaussiana en la energia
  Y distribucion Uniforme en el angulo */

  //Marsaglia algoritm
  rowvec x(dimension);
  double radio, angulo;
  double centro; 
  double delta;
  centro=sqrt(hbarra * nivel);
  delta=sqrt(hbarra * (nivel+1))-sqrt(hbarra * nivel);

  x.randn();
  x=x/norm(x,2);
  // aqui ya tenemos un punto en la esfera de radio 1

  radio=as_scalar(randu(1))*delta;
  radio=pow(radio+centro, 1.0/(double)dimension);
  radio=radio*pow(delta+centro, (double)(dimension-1)/(double)dimension);
  
  x=x*radio;
  
  return x;
  
};


mat Populate4BallShell(int nivel=0, int muestreo=1){
  //Crea un array del tamano de muestreo
  // poblando gausianamente cerca de la capa de energia que
  // nos interesa, dada por nivel.
  
  mat points(muestreo, 4);
  

  for(int i=0; i<muestreo; i++) {
    points.row(i)= CreateRandomPointinBall(nivel, 4);
  };

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



gsl_complex  WeylFunction2D
(int muestreo, simplectic * x, simplectic * y, 
 simplectic  mu, simplectic xhi){
  gsl_complex Weyl;
  
  double weylreal=0.00;
  double weylimag=0.00;

  for(int k=0; k<muestreo; k++){
    //Aqui hacemos la transformada de Fourier
	    // Que nos da la Funcion de Weyl
    
    weylreal+=cos((x[k].simplecticproduct(mu)+
		   y[k].simplecticproduct(xhi))/hbarra);
    
    weylimag+=sin((x[k].simplecticproduct(mu)+
		   y[k].simplecticproduct(xhi))/hbarra);
    
  };
  
  
  
  weylreal=weylreal/(hbarra*pi*muestreo);
  weylimag=weylimag/(hbarra*pi*muestreo);	  
  Weyl=gsl_complex_rect(weylreal, weylimag);
  
  cout<<mu.q<<"\t"<<mu.p<<"\t"
      <<xhi.q<<"\t"<<xhi.p<<"\t"
      <<weylreal<<"\t"<<weylimag
      <<endl;

  return Weyl;

};
