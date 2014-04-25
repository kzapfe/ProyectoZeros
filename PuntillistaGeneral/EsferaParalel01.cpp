//Probemos hacer gato estados cuadridimensionales.

#include <iostream>
#include <fstream>
#include <vector>
#include <armadillo>

#include <ctime> //para ver cuanto se tarda en cada punto
#include <omp.h> //para ver cuanto se tarda en cada punto

//Ultra Global Variables
const int maximumgauss=1400;
const int resol=100;
double magnituddeinteres=2.7;

#include "QuantumConstants.hpp"
#include "PopulateSphericalShells01.hpp"
#include "SeccionesPointillistas01.hpp"



using namespace std;
using namespace arma;

/********************************************************************/

int main(){
  //No mames que bien trabajas despues de hacer ejercicio 
  // Y SIN INTERNET. Recuerda eso SIEMPRE.

  
  //El numero de estados coherentes. El numero de gatos es 
  // la combinacion de pares posibles.

  //Los centros, otra vez para variar.
  mat centros;
  centros=Populate4BallShell(50, maximumgauss);
  int maximumgauss=centros.n_rows;
  //iniocilizar el estado poblado
  //parece ser que vector no funciona aqui bien

  vector   <simplectic> x,y;
  x.resize(maximumgauss);
  y.resize(maximumgauss);


  for(int i=0; i<maximumgauss; i++){
    //Selecciona centros en los subespacios X y Y
    x.at(i).SetQandP(centros(i,0), centros(i,1));
    y.at(i).SetQandP(centros(i,2), centros(i,3));
  };

  
  /*En esta aproximacion no necesitamos calcular
    wigner: es solamente dibujar con Gnuplot
    las diversas secciones del siguiente archivo */

  centros.save("CentroEsferas01.dat", arma_ascii);
  
  clock_t start, finish;
  start=clock();
  
#pragma omp parallel sections num_threads(4)
  {

      
  /*Tabla Booleana */
  //No te interesan dos cortes: los cruzados
  /* 0 0 -> qx qy 
     0 1 -> qx px
     1 0 -> qy py
     1 1 -> px py */
     




#pragma omp section 
    {
      bool qx,qy;
      string testeando;
      
      qx=false;
      qy=false;
      testeando="CorteExactDiracCuerdasQGrande.dat";

      WeylSection(x, y, maximumgauss,
		  qx, qy, testeando);		          

    }
     

#pragma omp section
    {
    
           bool qx,qy;
	   string testeando;
      
	   qx=false;
	   qy=true;
	   testeando="CorteExactDiracCuerdasMuGrande.dat";
	   WeylSection(x, y, maximumgauss,
		       qx, qy, testeando);		          
	   


    }
    
#pragma omp section
    {

      
           bool qx,qy;
	   string testeando;
      
	   qx=true;
	   qy=false;
	   testeando="CorteExactDiracCuerdasXhiGrande.dat";
	   WeylSection(x, y, maximumgauss,
			   qx, qy, testeando);		          
	       

    }

#pragma omp section
    {

      
           bool qx,qy;
	   string testeando;
      
	   qx=true;
	   qy=true;
	   testeando="CorteExactDiracCuerdasPGrande.dat";

	   WeylSection(x, y, maximumgauss,
		       qx, qy, testeando);		          
	       
   
      
    } //Cierra el ultimo  pragma omp section de cuerdas

  } //Cierra el pragma paralel omp sections directive

  finish=clock();


  cout<<"Nos tardamos "<<finish-start<< " Ciclos de cpu "<<endl; 
  
  return 0;
 
}

