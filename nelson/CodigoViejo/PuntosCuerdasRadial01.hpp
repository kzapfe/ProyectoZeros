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
  //En esta version mu esta repartido
  // en el plano rho theta, en lugar de p q
  
  //Solo piensa en superposicion de ondas medio aleatoreas
  const int resolucion=200; //La resolucion es una propiedad de la rutina.
  
  //Deberias poner esto en una rutina externa
  //Es una patada en los huevos

  std::ostringstream escupechord;
  std::ostringstream escupecondini;
  std::ostringstream escuperadial;
  std::ostringstream escupezeros;

  escupechord<<muestreo<<"_ChordsExtraDeform.dat"<<std::ends;
  escupecondini<<muestreo<<"_condinixtraRadial.dat"<<std::ends;
  escuperadial<<muestreo<<"_ChordsRadialHomogen.dat"<<std::ends;
  escupezeros<<muestreo<<"_Zeros.dat"<<std::ends;

  std::string stringchord;
  std::string stringcondini;
  std::string stringradial;
  std::string stringzeros;

  stringchord=escupechord.str();
  stringcondini=escupecondini.str();
  stringradial=escuperadial.str();
  stringzeros=escupezeros.str();

  const char *nombrechords=stringchord.c_str();
  const char *nombrecondini=stringcondini.c_str();
  const char *nombreradial=stringradial.c_str();
  const char *nombrezeros=stringzeros.c_str();


  ofstream chordspace;
  ofstream condinispace;
  ofstream chordradialspace;
  ofstream zeroline;


  
  chordspace.open(nombrechords);
  condinispace.open(nombrecondini);
  chordradialspace.open(nombreradial);
  zeroline.open(nombrezeros);



  
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
    	double OldWeylReal=0.0000000;
	double NewWeylReal=0.0000000;

    for(int j=0;j<(2*resolucion);j++){

      
	double extrem=0.6;
	double theta, rho;


	theta=pi*(double)i/(double)resolucion;
	rho=extrem*(double)j/(double)resolucion;

	mu.q= cos(theta)*rho;
	mu.p= sin(theta)*rho;
	  
	

	double weylreal=0.00;
	double weylimag=0.00;
	gsl_complex Weyl;
  
	

	for(int k=0; k<muestreo; k++){
	  weylreal+=cos((x[k].simplecticproduct(mu))/hbarra);
	  weylimag+=sin((x[k].simplecticproduct(mu))/hbarra);

	};

	weylreal=weylreal/(hbarra*pi*muestreo);
	weylimag=weylimag/(hbarra*pi*muestreo);
	NewWeylReal=weylreal;
	Weyl=gsl_complex_rect(weylreal, weylimag);

	chordspace<<mu.p<<"\t"<<mu.q<<"\t"
		  <<weylreal<<"\t"<<weylimag<<"\t"
		  <<gsl_complex_abs(Weyl)<<"\t"<<gsl_complex_arg(Weyl)
		  <<endl;
	gsl_complex zhi;
	 zhi=gsl_complex_rect(mu.p,mu.q);
	 chordradialspace<<rho<<"\t"
			 <<theta<<"\t"
			 <<weylreal<<"\t"<<weylimag<<"\t"
			 <<gsl_complex_abs(Weyl)<<"\t"
			 <<gsl_complex_arg(Weyl)<<endl;


	 //Aqui buscamos los Zeros

	 if(((NewWeylReal>0.00000000)&&(OldWeylReal<0.00000000))||
	    ((NewWeylReal<0.00000000)&&(OldWeylReal>0.00000000))){
	   zeroline<<rho<<"\t"
		   <<theta<<"\t"
		   <<weylreal<<"\t"<<weylimag<<"\t"
		   <<gsl_complex_abs(Weyl)<<"\t"
		   <<gsl_complex_arg(Weyl)<<endl;
	 };
	   

	 OldWeylReal=NewWeylReal;

      };
            chordspace<<endl;  
	    chordradialspace<<endl;  
  };


  
  


  chordspace.close();
  condinispace.close();
  chordradialspace.close();
  zeroline.close();

  return;

};
