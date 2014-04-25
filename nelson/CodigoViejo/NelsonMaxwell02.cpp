/*
Este programa dibuja la funcion de Cuerdas para 
planos mu=cte, en subdominio xhi. 

Version Maxwell 01: en RutinasNelson03 metemos una distribucion
que se sale de la camada de energia.
Para ello usaremos una rutina gsl 

Version 02:  vamos a meterle unos cuantos Deltas
Negativos, un cuarto aprox. 

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
  srand(893);

  //Parametros Globales
 
  cout<<"la seccion tiene q y p: "<<endl;
  cout<<mu.q<<"\t"<<mu.p<<endl;
  //Abrimos los archivos necesarios

  std::ostringstream escupeweyl;
  std::string hazstring;

  escupeweyl<<"NegTest"<<"_"<<
    "Chords_0-0.dat"<<std::ends;
  hazstring=escupeweyl.str();
  const char *nombreweyl=hazstring.c_str();
  ofstream WeylSeccion;
  WeylSeccion.open(nombreweyl);

 
  std::ostringstream escupemedias;

  escupemedias<<"NegTest"<<
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
  DiracDeltas=PopulateNelsonwithGauss(nivel, muestreo);
  
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


  


  //----------------------------------------------------
  
  //********HERE BEGINS THE PARTY !!!******************
  //----------------------------------------------------


      
  double argument;
  double azarneg;
  
  for(int i=-resolucion;i<resolucion;i++){
    for(int k=-resolucion;k<resolucion;k++){
      
	 
      double weylreal=0.00;
      double weylimag=0.00;
      simplectic xhi(0.0,0.0);
      
	     
	
      xhi.q=(double)i*extrem/(double)resolucion;
      xhi.p=(double)k*extrem/(double)resolucion;

      
	
	for(int k=0; k<muestreo; k++){
	  //Aqui hacemos la transformada de Fourier
	  // Que nos da la Funcion de Weyl
	  
	  azarneg=as_scalar(randu(1));
	  // cout<<" eeee azar negativo eeee "<<azarneg<<endl;

	  argument=-(x[k].simplecticproduct(mu)+
		     y[k].simplecticproduct(xhi)-4.0*xhi.p)/hbarra;
	  

	  if(azarneg<0.25){
	    //Negative Deltas
	    	  weylreal-=cos(argument);	  
		  weylimag-=sin(argument);
	  }else{
	    //Positive Deltas
	  weylreal+=cos(argument);	  
	  weylimag+=sin(argument);
	  };
	  
	};
	
	     
       	weylreal=weylreal/sqrt(hbarra*pi*muestreo);
	weylimag=weylimag/sqrt(hbarra*pi*muestreo);
	
	WeylSeccion<<xhi.q<<"\t"<<xhi.p<<
	  "\t"<<weylreal<<"\t"<<weylimag<<endl;
	
		
	
    };
    
    
    WeylSeccion<<endl;
    
    
    cout<<"Vamos en la linea i="<<i<<endl;
  
  };

  
  WeylSeccion.close();
  
 
  return 0;
  
}
