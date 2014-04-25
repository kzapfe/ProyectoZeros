/*
Programa que obtiene la funcion de Husimi a partir 


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

  //Matriz de Wigner from Fabrizio
  //La original tiene lineas blancas, aqui la tienes que poner sin ellas.
  mat DatWigner;
  DatWigner.load("DatosLimpios.dat",raw_ascii);

  //elegante seria cargar esto del otro programa
  // o de los datos, pero hoy tienes prisa
  //nos choteamos los valores del otro programa 
  double deltap,deltaq;
  deltaq=0.20/400.0;
  deltap=0.20/400.0;
  
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

  cout<<"deltaq="<<deltaq<<"\t"<<"deltap="<<deltap<<"\n";

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
	Husimi+= DatWigner(i,2)
	  //exp(-1/hbarra)
	  * exp((-pow((DatWigner(i,0)-Equis.q),2)-pow((DatWigner(i,1)-Equis.p),2))/hbarra)
	  *deltap*deltaq;
	/*cout<<"Llevamos esto de Husimi "
	    <<Equis.q<<"\t"
		    <<Equis.p<<"\t"
		    <<Husimi<<endl;   */ 
	
      };
  
      Husimi=Husimi/(pi*hbarra);

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
