//Probemos hacer gato estados cuadridimensionales.

#include <iostream>
#include <fstream>
#include <vector>
#include <armadillo>

#include <ctime> //para ver cuanto se tarda en cada punto
#include <omp.h> //para ver cuanto se tarda en cada punto

//Ultra Global Variables


#include "simplectic01.hpp"
#include "RutinaExternaEstupida01.hpp"



using namespace std;
using namespace arma;

/********************************************************************/

int main(int argc, char* argv[]){
  //No mames que bien trabajas despues de hacer ejercicio 
  // Y SIN INTERNET. Recuerda eso SIEMPRE.

  string NombreCentros;
  NombreCentros=argv[1];
  //Los centros, otra vez para variar.
  mat centros;
  centros.load(NombreCentros);
  int maximumgauss=centros.n_rows;

  //iniocilizar el estado poblado
  //parece ser que vector no funciona aqui bien

  vector   <simplectic> x,y;
  simplectic auxX, auxY;
  x.resize(maximumgauss);
  y.resize(maximumgauss);

  for(int i=0; i<maximumgauss; i++){
    //Selecciona centros en los subespacios X y Y
 
    x.at(i).SetQandP(centros(i,0), centros(i,1));
    y.at(i).SetQandP(centros(i,2), centros(i,3));
    
    //cout<<x[i].q<<endl;
  };
  
  cout<<"seg fault ? "<<"juat"<< endl;

    
  string testeando;
  testeando="PasameUnVector.dat";
  
  TestRecibeVector(x, y,testeando);		          

 
  
  return 0;
 
}

