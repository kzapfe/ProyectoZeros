//Veamos que podemos hacer con esto
//Tratemos de encontrar los Blind Spots Subplanckianos
//con esta chiva.


//No necesitas incluir nada que ya se haya incluido antes
#include "RutinasExternas01.hpp"
#include <fstream>
#include <iostream>




using namespace std;

//Esto es pleonasmo, pero para que te acuerdes.
using namespace arma;


int main(){

  srand(23890);

  //Solo piensa en superposicion de ondas medio aleatoreas
  const int muestreo=500;
  const int nivel=300;
  const int resolucion=200;
  ofstream chordspace;
  
  chordspace.open("Cuerdaspacezero03.dat");

  
  //Bueno, si el sistema es caotico 
  //tiene al menos dos grados de libertad
  simplectic * x;
  simplectic * y;
  simplectic xhi(0.1,0.1), mu(0.000,0.000);

  
  x=new simplectic[muestreo];
  // y=new simplectic[1];
  
 
  x=PopulateEnergyLevel( nivel, muestreo);
  
  
  
  for(int i=-resolucion;i<resolucion;i++){
      for(int j=-resolucion;j<resolucion;j++){
	
	//nos interesa hasta el primer Zero
	double extrem=0.3;
	
	mu.q= extrem *
	  (double)i/(double)resolucion;
	mu.p= extrem *
	  (double)j/(double)resolucion;


	double weylreal=0.00;
	double weylimag=0.00;
  
	for(int k=0; k<muestreo; k++){

	  weylreal+=cos((x[k].simplecticproduct(mu))/hbarra);
	  weylimag+=sin((x[k].simplecticproduct(mu))/hbarra);

	};

	weylreal=weylreal/(hbarra*pi*muestreo);
	weylimag=weylimag/(hbarra*pi*muestreo);

  
	chordspace<<mu.p<<"\t"<<mu.q<<"\t"
		  <<weylreal<<"\t"<<weylimag<<endl;
      };
      chordspace<<endl;  
  };


  chordspace.close();


}
