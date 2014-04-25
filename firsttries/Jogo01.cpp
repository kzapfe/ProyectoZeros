//Veamos que podemos hacer con esto



#include <fstream>
#include <cstdlib>
#include <iostream>
#include "simplectic01.hpp"
#include "coherentstate01.hpp"

using namespace std;


int main(){
  
  ofstream wigner, cordas;
  
  simplectic chi(1.0,0.0);
  coherent estado1(chi);
  simplectic eta(cos(pi*2.0/3.),sin(pi*2.0/3.0));
  coherent estado2(eta);
  simplectic mu(cos(pi*4.0/3.),sin(pi*4.0/3.0));
  coherent estado3(mu);

  cout<<eta.q<<"\t"<<eta.p<<endl;
  cout<<mu.q<<"\t"<<mu.p<<endl;


  cout<<endl;
  cout<<estado1.zhi.q<<"\t"<<estado1.zhi.p<<endl;  
  cout<<estado2.zhi.q<<"\t"<<estado2.zhi.p<<endl;  
  cout<<estado3.zhi.q<<"\t"<<estado3.zhi.p<<endl;


  

double varx,vary;

  

  wigner.open("Wigner.dat");
  cordas.open("cordas.dat");
  
  //cout<<chi.p<<"\t"<<chi.q<<endl;
  
  
  
  

  for(int i=-150;i<150;i++){
    varx=(double)i/75.0;
    gsl_complex amplitud1,  amplitud2;
    /*    
    amplitud1=estado1.qspace(varx);
    amplitud2=estado1.pspace(varx);
    
    cout<<varx<<"\t"<<
      GSL_REAL(amplitud1)<<"\t"<<
      GSL_IMAG(amplitud1)<<"\t"<<
      GSL_REAL(amplitud2)<<"\t"<<
      GSL_IMAG(amplitud2)<<endl;
    */

    for (int j=-150;j<150;j++){
      vary=(double)j/75.0;
      simplectic y(varx,vary);
      double wig;
      wig=estado1.Wigner(y)+estado2.Wigner(y)+estado3.Wigner(y);
      
      gsl_complex Xi;

      Xi=gsl_complex_add(estado1.cuerdas(y),estado2.cuerdas(y));
      Xi=gsl_complex_add(Xi,estado3.cuerdas(y));

      wigner<<varx<<"\t"<<
	vary<<"\t"<<wig<<endl;

      cordas<<varx<<"\t"<<
	vary<<"\t"<<GSL_REAL(Xi)
	    <<"\t"<<GSL_IMAG(Xi)
	    <<endl;

    
    };

    
        
  };
  
  wigner.close();
  cordas.close();

  return 0;

}
