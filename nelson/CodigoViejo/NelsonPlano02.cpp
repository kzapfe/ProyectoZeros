/*
Este programa dibuja la funcion de Cuerdas para 
planos mu=cte, en subdominio xhi. 

Version 02: Vamos a poner este desmadre en orden
 */
#include <armadillo>
#include "simplectic01.hpp"
#include "ParametrosGlobales.hpp"
#include "RutinasNelson02.hpp"


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

  escupeweyl<<"HiRes"<<"_"<<"08138_"<<
    "Chords_0-0.dat"<<std::ends;
  hazstring=escupeweyl.str();
  const char *nombreweyl=hazstring.c_str();
  ofstream WeylSeccion;
  WeylSeccion.open(nombreweyl);

 

  std::ostringstream escupemedias;

  escupemedias<<"HiRes"<<"_"<<"08138"<<
    "_MediasWigner0-0.dat"<<std::ends;
  hazstring=escupemedias.str();
  const char *nombremedwig=hazstring.c_str();
  ofstream MediasWigner;
  MediasWigner.open(nombremedwig);

  ofstream Centros;
  Centros.open("CentrosDiracOtraVez.dat");



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
  

  for(int i=0; i<muestreo/10;i++){
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
	  //La convencion que estas usando es el producto ^ a la Alfredo
	  //x.simplecticproduct(mu)= x^mu=x.p*mu.q-x.q*mu.p
	  // y entonces la FT que estas usando es
	  // Integral( Wigner (x) *exp (-i/h *x^mu))

	  weylreal+=cos(-(x[k].simplecticproduct(mu)+
			 y[k].simplecticproduct(xhi))/hbarra);
	  
	  weylimag+=sin(-(x[k].simplecticproduct(mu)+
			  y[k].simplecticproduct(xhi))/hbarra);
	  
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
