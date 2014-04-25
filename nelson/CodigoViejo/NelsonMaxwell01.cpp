/*
Este programa dibuja la funcion de Cuerdas para 
planos mu=cte, en subdominio xhi. 

Version Maxwell 01: en RutinasNelson03 metemos una distribucion
que se sale de la camada de energia.
Para ello usaremos una rutina gsl 


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

//Esto es pleonasmo, pero para que te acuerdes.

int main(){
  
  //inicializar el semillador
  srand(8931);

  //Parametros Globales
 
  cout<<"la seccion tiene q y p: "<<endl;
  cout<<mu.q<<"\t"<<mu.p<<endl;
  //Abrimos los archivos necesarios

  std::ostringstream escupeweyl;
  std::string hazstring;

  escupeweyl<<"Xhi_P_Qequals-0-0"<<"_"<<"08138"<<
    "_Chords.dat"<<std::ends;
  hazstring=escupeweyl.str();
  const char *nombreweyl=hazstring.c_str();
  ofstream WeylSeccion;
  WeylSeccion.open(nombreweyl);

 
  std::ostringstream escupemedias;

  escupemedias<<"CorrectFourierK"<<"_"<<"08138"<<
    "_MediasWigner0-0.dat"<<std::ends;
  hazstring=escupemedias.str();
  const char *nombremedwig=hazstring.c_str();
  ofstream MediasWigner;
  MediasWigner.open(nombremedwig);

  ofstream Centros;
  Centros.open("ExactNCentrosDirac.dat");


  //----------------------------------------------------
  //Poblar la camada de energia
  mat DiracDeltas;
  DiracDeltas=PopulateNelson(nivel, muestreo);
  
  simplectic * x;
  simplectic * y;
  
  x=new simplectic[muestreo];
  y=new simplectic[muestreo];


  for(int i=0;i<muestreo;i++){
    x[i].q=DiracDeltas(i,0);
    x[i].p=DiracDeltas(i,1);
    y[i].q=DiracDeltas(i,2);
    y[i].p=DiracDeltas(i,3);
  };
  

  for(int i=0; i<muestreo;i++){
    //Cheq that everything is as it should
    Centros<<x[i].q<<"\t"<<x[i].p<<"\t"
	   <<y[i].q<<"\t"<<y[i].p<<"\t"
		<<NelsonEnergy(DiracDeltas.row(i))<<endl;
  };


  cout<<" Acabamos de popular la capa de energia"
      << " con "<<muestreo<<" puntos"<<endl;
  
  mat meanvalues;
  
  meanvalues=mean(DiracDeltas);
  
  MediasWigner<<"Primer momento x_q: "<< meanvalues(0)<<endl;
  MediasWigner<<"Primer momento x_p: "<< meanvalues(1)<<endl;
  MediasWigner<<"Primer momento y_q: "<< meanvalues(2)<<endl;
  MediasWigner<<"Primer momento y_p: "<< meanvalues(3)<<endl;
  MediasWigner.close();

  Centros.close();
  


  //----------------------------------------------------
  
  //********HERE BEGINS THE PARTY !!!******************
  //----------------------------------------------------


      
  double argument;
  double intinteg;
  //La longitud del intervalo de integracion
  intinteg=extrem/(double)resolucion;

  for(int i=-resolucion;i<resolucion;i++){
    for(int k=-resolucion;k<resolucion;k++){
      	 
      double weylreal=0.00;
      double weylimag=0.00;
      simplectic xhi(0.0,0.0);
      
      //Cambiamos la secccion a las cuerdas_p
      xhi.q=0.0;
      mu.p=(double)i*intinteg;
      xhi.p=(double)k*intinteg;
	
  	for(int k=0; k<muestreo; k++){
  	  //Aqui hacemos la transformada de Fourier
  	  // Que nos da la Funcion de Weyl
	  
  	  argument=+(x[k].simplecticproduct(mu)+
  		     y[k].simplecticproduct(xhi))
  	    /hbarra;
	  
  	  weylreal+=cos(argument);
	  
  	  weylimag+=sin(argument);
	  
  	};
	
  	//esta constante de normalizacion es incorrecta.
  	// la correcta es la constante "1/muestreo"
       	weylreal=weylreal/sqrt(hbarra*pi)*intinteg*intinteg;
  	weylimag=weylimag/sqrt(hbarra*pi)*intinteg*intinteg;
	
  	WeylSeccion<<xhi.p<<"\t"<<mu.p<<
  	  "\t"<<weylreal<<"\t"<<weylimag<<endl;
			
	
    };
    
    
    WeylSeccion<<endl;
    
    
    cout<<"Vamos en la linea i="<<i<<endl;
  
  };

  
  WeylSeccion.close();
  
 
  return 0;
  
}
