//Veamos que podemos hacer con esto
//Tratemos de encontrar los Blind Spots Subplanckianos
//con esta chiva.


//No necesitas incluir nada que ya se haya incluido antes
#include "RutinasExternas01.hpp"
#include <fstream>
#include <iostream>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>




using namespace std;

//Esto es pleonasmo, pero para que te acuerdes.
using namespace arma;


void CalculaWeyl(int muestreo, int nivel ){

  
  //Solo piensa en superposicion de ondas medio aleatoreas
  const int resolucion=200; //La resolucion es una propiedad de la rutina.
  

  std::ostringstream escupechord;
  std::ostringstream escupecondini;

  escupechord<<muestreo<<"_Chords.dat"<<std::ends;
  escupecondini<<muestreo<<"_condini.dat"<<std::ends;

  std::string stringchord;
  std::string stringcondini;
  stringchord=escupechord.str();
  stringcondini=escupecondini.str();

  const char *nombrechords=stringchord.c_str();
  const char *nombrecondini=stringcondini.c_str();


  ofstream chordspace;
  ofstream condinispace;
  
  chordspace.open(nombrechords);
  condinispace.open(nombrecondini);




  
  //Bueno, si el sistema es caotico 
  //tiene al menos dos grados de libertad
  simplectic * x;
  simplectic * y;
  simplectic xhi(0.1,0.1), mu(0.000,0.000);

  
  x=new simplectic[muestreo];
  // y=new simplectic[1];
  
 
  x=PopulateEnergyLevel( nivel, muestreo);
  
  for(int i=0;i<2756;i++){
    condinispace<<x[i].q<<"\t"<<x[i].p<<endl;
      };
  
  for(int i=-resolucion;i<resolucion;i++){
      for(int j=-resolucion;j<resolucion;j++){
	double extrem=0.6;
	
	mu.q= extrem *
	  (double)i/(double)resolucion;
	mu.p= extrem *
	  (double)j/(double)resolucion;


	double weylreal=0.00;
	double weylimag=0.00;
	gsl_complex Weyl;
  
	

	for(int k=0; k<muestreo; k++){

	  weylreal+=cos((x[k].simplecticproduct(mu))/hbarra);
	  weylimag+=sin((x[k].simplecticproduct(mu))/hbarra);

	};

	weylreal=weylreal/(hbarra*pi*muestreo);
	weylimag=weylimag/(hbarra*pi*muestreo);

	Weyl=gsl_complex_rect(weylreal, weylimag);

	chordspace<<mu.p<<"\t"<<mu.q<<"\t"
		  <<weylreal<<"\t"<<weylimag<<"\t"
		  <<gsl_complex_abs(Weyl)<<"\t"<<gsl_complex_arg(Weyl)
		  <<endl;
      };
      chordspace<<endl;  
  };


  chordspace.close();
  condinispace.close();

  return;

};
