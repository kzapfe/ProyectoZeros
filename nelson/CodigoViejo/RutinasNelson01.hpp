#include <armadillo>
#include "simplectic01.hpp"


using namespace arma;

//Constantes del programa
// y del Universo.
double const pi=3.14151926535898; 
double const hbarra=1.000;



double NelsonPotential(rowvec x){

  double result;
  double q1,q2,p1,p2;
  
  q1=x(0);
  p1=x(1);
  q2=x(2);
  p2=x(3);
  
   result=(p1*p1+p2*p2)/2.0+
    q1*q1/20.0
    +(q2-q1*q1/2)*(q2-q1*q1/2);
   
  return result;

};


double ValorExactoq2(double q1, double  Potencial){
  
  //Aqui tambien tiramos un dado
  double result;
  double aux=0.000;
  aux=as_scalar(randn(1));
  if(aux>0.5000){
    result=(q1*q1+sqrt(4.0*Potencial-q1*q1/20.0))/2.0;
  }else{
    result=(q1*q1-sqrt(4.0*Potencial-q1*q1/20.0))/2.0;
  };
  
  return result;

};

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


mat PopulateNelson(int nivel=0, int muestreo=1){
  //Usa la cabezota, guey
  //para cada valor de E, q2 y |p| pueden ser exactamente determinados.
  //entonces realmente lo unico que necesitamos
  //es escoger aleatoreamente TRES de las cuatro, y variarle una
  //unidad a E.

  double particion;
  double Energia, Potencial, Kinetica;
  mat result;
  double Ham;
  int i=0;
  

  //Tu unidad de energia mas pequenha es 0.05*hbarra (freq. menor)
  //First Coin Toss: lugar en el gordito entre nivel y  nivel+1 .
  Energia=(double)nivel+as_scalar(randu(1))*hbarra*0.05;
  

  //Second Coin toss: cuanto le toca a Kin y cuanto a Pot
  
  
  result=zeros<mat>(muestreo, 4);
  

  for(int i=0; i<muestreo;i++){
    
    Kinetica=Energia*as_scalar(randu(1));
    Potencial=Energia-Kinetica;
  

    rowvec auxiliar(4);
    double theta;
    //Third Coin Toss: direccion del momento
    theta=2.00*pi*as_scalar(randu(1));
        
    
    //momentos 1 y 2
    auxiliar(1)=sqrt(2.0*Kinetica)*cos(theta);
    
    auxiliar(3)=sqrt(2.0*Kinetica)*sin(theta);
    

    //esto son posiciones


    auxiliar(0)=as_scalar(randu(1))*
      (2.0*sqrt(20.0*Potencial))-sqrt(20.0*Potencial);
    auxiliar(2)=ValorExactoq2(auxiliar(0), Potencial);

    result.row(i)=auxiliar;
   
  };
      
     

  return result;

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


