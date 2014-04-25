
#include "simplectic01.hpp"
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <cmath>

#include "coherentstate01.hpp"



# ifndef __coherentsum__
# define __coherentsum__

class coherentsum{

private:
  //Por alguna puta razon no funciona mas que en private...
  int totals;
  gsl_complex pesos[];
  simplectic centros[]; 

public:

  
  coherent componente[];


  coherentsum(simplectic center[], gsl_complex wheights[], int &chan)
  {
    /*Esto construye simplemente de forma abstracta la suma de estados
      coherentes*/
   
    /*Tenemos que dar el numero de componentes desde antes, por la pendejada
     del c y c++ de no distinguir entre bla[]
     */
    totals=*&chan;
    

   
    for(int i=0; i<totals; ++i){
    
      centros[i] = center[i];  
      pesos[i] = wheights[i];
      
      //demasiado rollo pero parece funcionar
      simplectic  aux;
      aux.q=center[i].q;
      aux.q=center[i].p;
      coherent aux2(aux);
      componente[i]=aux2;
      cout<<GSL_REAL(pesos[i])<<endl;
    };
 

    

  };//Here it ends the constructor

  




///-------------------------------666---------------------------
  /// SHUT UP AND TAKE MY MONEY!!!----
  /// Here shall began the methods!
  ///-------------------------------666---------------------------




  gsl_complex qspace(double x=0.00){
    //La funcion de onda en el espacio q
    // x aqui es una  posicion
    //heredamos todo de coherentstate
    gsl_complex result;
    result=gsl_complex_rect(0.00,0.000);
    for(int i=0; i<totals; i++){
      result=
	gsl_complex_add(result, 
			gsl_complex_mul(componente[i].qspace(x),pesos[i]));
      
    };

    return result;

  };


 gsl_complex pspace(double x=0.00){
    //La funcion de onda en el espacio p
    // x aqui es un momento
    //heredamos todo de coherentstate
    gsl_complex result;
    result=gsl_complex_rect(0.00,0.000);
    
    for(int i=0; i<totals; i++){
      result=
	gsl_complex_add(result, 
			gsl_complex_mul(componente[i].pspace(x),pesos[i]));
    };

    return result;

  };






};


#endif
