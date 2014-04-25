/*
Programa que obtiene la funcion de Husimi a partir 
de los datos de Karel. Como aqui la funcion de Wigner
siempre tiene el valor 1 donde es distinta de cero,
el codigo es mas simple que en el caso de Fabricio,
ya que solo integramos Deltas de Dirac
 */
// estos dos son los de gsl para rutinas aleatoreas
//#include <gsl_rng.h>
//#include <gsl_randist.h>
// De momento usaremos solo la guassiana de armadillo

#include <cstdlib>
#include <armadillo>
#include "simplectic01.hpp"
#include "ParametrosGlobales.hpp"



using namespace std;
using namespace arma;

//Esto es pleonasmo, pero para que te acuerdes.

int main(){

  //Matriz de Wigner from Karel
  //La original tiene 4 columnas, estamos interesados en la 3 y 4, notacion armadillos
  // son la 2 y 3
  mat DatWigner;
  DatWigner.load("CentrosWigner.dat",raw_ascii);

    
  
  const double maxequisp=1.500;
  const double minequisq=-1.500;
  const double maxequisq=4.00;
  const int resolucion=200;

  double Husimi;
  simplectic Equis;

  ofstream HusimiSection;
  HusimiSection.open("HusimiKarel01.dat");

  int n_datos;
  n_datos=DatWigner.n_rows;

 

  for(int l=0.00; l<=2.0*resolucion; l++){ 
    //aqui arriba estamos haciendo que las l tengan la
    //misma res que las k.
    
    for(int k=-resolucion; k<=resolucion; k++){ 
      
      Husimi=0.00;
      
      //esta parte requiere otra forma de repartir bien los datos
      Equis.q=minequisq+
	(double)l/((double)resolucion*2)*
	(maxequisq-minequisq);
      
      //esta parte esta fina asi, va de negativo a positivo
      Equis.p=(double)k/((double)resolucion)*maxequisp;

     
      
      // Esta parte es la convolucion con los datos de la funcion de Wigner.
      for(int i=0; i<n_datos; i++){
	//El lugar 2 y 3 de mis datos son y_q y y_p del Centro de Dirac. 
	Husimi+=
	  exp((-pow((DatWigner(i,2)-Equis.q),2)-pow((DatWigner(i,3)-Equis.p),2))/hbarra);
	 	
      };


      Husimi=Husimi/(pi*hbarra*n_datos);

      HusimiSection<<Equis.q<<"\t"
		   <<Equis.p<<"\t"
		   <<Husimi<<endl;    
      
      /*
      cout<<Equis.q<<"\t"
		    <<Equis.p<<"\t"
		    <<Husimi<<endl;    
      
      */
	    
    };

    cout<<"vamos en el l="<<l<< " de lmax="<<2 *resolucion<<endl;
  };

}
