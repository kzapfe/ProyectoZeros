/*
Este programa dibuja la funcion de Cuerdas para 
planos mu=cte, en subdominio xhi. 

Version 02: Vamos a poner este desmadre en orden
 */
#include <armadillo>
#include "simplectic01.hpp"
#include "ParametrosGlobales.hpp"


using namespace std;
using namespace arma;
//Esto es pleonasmo, pero para que te acuerdes.

  
int main(int argc, char *argv[]){
 
  cout<<"la seccion tiene q y p: "<<endl;
  cout<<mu.q<<"\t"<<mu.p<<endl;
  //Abrimos los archivos necesarios

  ofstream WeylSeccion;
  WeylSeccion.open("Cuerdas.dat");

  ofstream MediasWigner;
  MediasWigner.open("Promedios.dat");

  //----------------------------------------------------
  //Poblar la camada de energia
  mat DiracDeltas;

  cout<<argv[1]<<endl;
  
  string centrosexterior=argv[1];
  
  //cout<<"ya la cagaste"<<endl;

  DiracDeltas.load(centrosexterior);
  
  simplectic * x;
  simplectic * y;

  int muestreo=DiracDeltas.n_rows;
  
  x=new simplectic[muestreo];
  y=new simplectic[muestreo];


  for(int i=0;i<muestreo;i++){
    x[i].q=DiracDeltas(i,0);
    x[i].p=DiracDeltas(i,1);
    y[i].q=DiracDeltas(i,2);
    y[i].p=DiracDeltas(i,3);
  };
  
  
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
	
	     
       	weylreal=weylreal/sqrt(hbarra*2.0*pi)/muestreo;
	weylimag=weylimag/sqrt(hbarra*2.0*pi)/muestreo;
	
	WeylSeccion<<xhi.q<<"\t"<<xhi.p<<
	  "\t"<<weylreal<<"\t"<<weylimag<<endl;
	
		
	
    };
    
    
    WeylSeccion<<endl;
    
    
    //cout<<"Vamos en la linea i="<<i<<endl;
  
  };

  
  WeylSeccion.close();
  
 
  return 0;
  
}
