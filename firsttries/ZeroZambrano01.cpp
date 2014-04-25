//Veamos que podemos hacer con esto
//Tratemos de encontrar los Blind Spots Subplanckianos
//con esta mierda.


#include <fstream>
#include <cstdlib>
#include <iostream>
#include "simplectic01.hpp"
#include "coherentstate01.hpp"

using namespace std;


int main(){
  


  ofstream correlation;
  int maxpuntos=200;
  
  double dispstep=0.15/(double)maxpuntos;



  //El estado inicial es el estado triangular simetrico
  simplectic chi(1.0,0.0);
  coherent estado1(chi);
  simplectic eta(cos(pi*2.0/3.),sin(pi*2.0/3.0));
  coherent estado2(eta);
  simplectic mu(cos(pi*4.0/3.),sin(pi*4.0/3.0));
  coherent estado3(mu);
  


  
  simplectic chi2(chi);
  simplectic eta2(eta);
  simplectic mu2(mu);
  

  double varx,vary;

  

  correlation.open("Correlation.dat");
  
  
  //cout<<chi.p<<"\t"<<chi.q<<endl;


  for(int k=-maxpuntos;k<maxpuntos;k++){
    for(int l=-maxpuntos;l<maxpuntos;l++){
  
      double vartransq,vartransp;
      vartransq=k*dispstep;
      vartransp=l*dispstep;
  
      chi2.q=vartransq+chi.q;
      eta2.q=vartransq+eta.q;
      mu2.q=vartransq+mu.q;

      chi2.p=vartransp+chi.p;
      eta2.p=vartransp+eta.p;
      mu2.p=vartransp+mu.p;
      

      coherent estadodesplazado1(chi2);   
      coherent estadodesplazado2(eta2);
      coherent estadodesplazado3(mu2);

      double wigcorr=0.0;



      // Este es el ciclo (doble) de las funciones
      //Aqui obtenemos el producto interno
      for(int i=-maxpuntos;i<maxpuntos;i++){
	double limiteinteres=2.0;
	
	varx=limiteinteres/(double)maxpuntos;
	
	for (int j=-maxpuntos;j<maxpuntos;j++){
	  vary=limiteinteres/(double)maxpuntos;
	    
	  simplectic y((i*varx),(j*vary));
	  double wignerzero,wignerdesplazada;
	  
	  wignerzero=estado1.Wigner(y)+estado2.Wigner(y)+estado3.Wigner(y);
	  wignerdesplazada=estadodesplazado1.Wigner(y)+
	    estadodesplazado2.Wigner(y)+estadodesplazado3.Wigner(y);
    
	  wigcorr+=wignerzero*wignerdesplazada*varx*vary;


	};//Cierra j 
      }; //Cierra i
      //Esto cierra el ciclo de integracion


      correlation<<vartransq<<"\t"<<
	vartransp<<"\t"<<wigcorr
	    <<endl;


    }; //cierra l

    cout<<"Llevo por ahi de "<<k<< " y voy hasta "<<maxpuntos<<endl;
  };//cierra k
  

  correlation.close();


  return 0;

}
