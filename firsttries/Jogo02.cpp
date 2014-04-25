//Veamos que podemos hacer con esto
//Intentemos hacer una superposicion
// de estados coherentes


#include <fstream>
#include <cstdlib>
#include <iostream>
#include "simplectic01.hpp"
#include "coherentstate01.hpp"
#include "coherentsuperposition.hpp"


using namespace std;


int main(){
  
  ofstream wigner, cordas;
  ofstream espacioq, espaciop;

  simplectic chi(1.0,0.0);
  simplectic eta(cos(pi*2.0/3.),sin(pi*2.0/3.0));
  simplectic mu(cos(pi*4.0/3.),sin(pi*4.0/3.0));
  
  int tantos=3;
  simplectic * centros=new simplectic[tantos];
  //gsl_complex * pesos=new gsl_complex[3];
  gsl_complex pesos[tantos];

  //este tipo de declaracion se llama array dinamico
  // por alguna razon el otro no FUNCIONA para nuestros propositos.
  
  
  
  centros[0]=chi;
  centros[1]=eta;
  centros[2]=mu;

  
  pesos[0]=gsl_complex_rect (1.0, 0.0);
  pesos[1]=gsl_complex_rect (-1.0, 0.0);
  pesos[2]=gsl_complex_rect (0.0, 1.0);

  
  coherentsum estado(centros, pesos, tantos);



  espacioq.open("EspacioPos.dat");
  espaciop.open("EspacioMom.dat");

  
  for(int i=-500 ; i<=500; i++){
    double aux;
    gsl_complex phi,psi;
    aux=(double)i/100.00;
  
    phi=estado.pspace(aux);
    psi=estado.qspace(aux);
    espacioq<<aux<<"\t"<<GSL_REAL(psi)<<"\t"<<GSL_IMAG(psi)<<endl;
    espaciop<<aux<<"\t"<<GSL_REAL(phi)<<"\t"<<GSL_IMAG(phi)<<endl;


  };
  

  //Esta ultima mierda es una cuestion de elegancia
  //  delete [] centros;
  //  delete [] pesos;

  cout <<"PIIIITO"<<endl;

  espaciop.close();
  espacioq.close();

  return 0;




 

}
