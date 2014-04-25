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

#include "CatStates01.hpp"
#include "PopulateSphericalShells01.hpp"
#include "CortesGatosBidimensionales01.hpp"


using namespace std;
using namespace arma;

/********************************************************************/

int main(){
  //No mames que bien trabajas despues de hacer ejercicio 
  // Y SIN INTERNET. Recuerda eso SIEMPRE.

  vector   <simplectic> x,y;
  //El numero de estados coherentes. El numero de gatos es 
  // la combinacion de pares posibles.

  //Los centros, otra vez para variar.
  mat centros;
  centros=Populate4BallShell(50, maximumgauss);


  //iniocilizar el estado poblado
  //parece ser que vector no funciona aqui bien
  CoherentState *CentroX, *CentroY;
  CatState *GatosX, *GatosY;

  CentroX=new CoherentState[maximumgauss];
  CentroY=new CoherentState[maximumgauss];

  GatosX=new CatState[maximumgauss*(maximumgauss-1)/2];
  GatosY=new CatState[maximumgauss*(maximumgauss-1)/2];



  for(int i=0; i<maximumgauss; i++){
    //Selecciona centros en los subespacios X y Y
    CentroX[i].SetCentre(centros(i,0), centros(i,1));
    CentroY[i].SetCentre(centros(i,2), centros(i,3));


  };

  int cuentainterferencias=0;

  for(int i=0; i<maximumgauss-1; i++){
    for(int j=i+1; j<maximumgauss; j++){
      //las interferencias son Indepèndientes
      GatosX[cuentainterferencias].SetGaussians(CentroX[i], CentroX[j]);
      GatosY[cuentainterferencias].SetGaussians(CentroY[i], CentroY[j]);      
      cuentainterferencias++;      
    }
  }


  //ahora toca calcular el valor de las chivas en la funcion:
  //Corte en el plano de las q y en el de las y.

  clock_t start, finish;
  start=clock();
  
#pragma omp parallel sections num_threads(8)
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
      testeando="CorteExactGaussQGrande.dat";

      WignerSection(CentroX, CentroY, GatosX, GatosY, 
		    maximumgauss, cuentainterferencias,
		    qx, qy, testeando);		          

    }
     

#pragma omp section
    {
    
           bool qx,qy;
	   string testeando;
      
	   qx=false;
	   qy=true;
	   testeando="CorteExactGaussXGrande.dat";

	   WignerSection(CentroX, CentroY, GatosX, GatosY, 
		    maximumgauss, cuentainterferencias,
		    qx, qy, testeando);


    }
    
#pragma omp section
    {

      
           bool qx,qy;
	   string testeando;
      
	   qx=true;
	   qy=false;
	   testeando="CorteExactGaussYGrande.dat";

	   WignerSection(CentroX, CentroY, GatosX, GatosY, 
		    maximumgauss, cuentainterferencias,
		    qx, qy, testeando);

    }

#pragma omp section
    {

      
           bool qx,qy;
	   string testeando;
      
	   qx=true;
	   qy=true;
	   testeando="CorteExactGaussPGrande.dat";

	   WignerSection(CentroX, CentroY, GatosX, GatosY, 
		    maximumgauss, cuentainterferencias,
		    qx, qy, testeando);

   
      
    } //Cierra el ultimo  pragma omp section de Wigner



#pragma omp section 
    {
      bool qx,qy;
      string testeando;
      
      qx=false;
      qy=false;
      testeando="CorteExactGaussCuerdasQGrande.dat";

      WeylSection(CentroX, CentroY, GatosX, GatosY, 
		    maximumgauss, cuentainterferencias,
		    qx, qy, testeando);		          

    }
     

#pragma omp section
    {
    
           bool qx,qy;
	   string testeando;
      
	   qx=false;
	   qy=true;
	   testeando="CorteExactGaussCuerdasMuGrande.dat";

	   WeylSection(CentroX, CentroY, GatosX, GatosY, 
		    maximumgauss, cuentainterferencias,
		    qx, qy, testeando);


    }
    
#pragma omp section
    {

      
           bool qx,qy;
	   string testeando;
      
	   qx=true;
	   qy=false;
	   testeando="CorteExactGaussCuerdasXhiGrande.dat";

	   WeylSection(CentroX, CentroY, GatosX, GatosY, 
		    maximumgauss, cuentainterferencias,
		    qx, qy, testeando);

    }

#pragma omp section
    {

      
           bool qx,qy;
	   string testeando;
      
	   qx=true;
	   qy=true;
	   testeando="CorteExactGaussCuerdasPGrande.dat";

	   WeylSection(CentroX, CentroY, GatosX, GatosY, 
		    maximumgauss, cuentainterferencias,
		    qx, qy, testeando);

   
      
    } //Cierra el ultimo  pragma omp section de cuerdas

  } //Cierra el pragma paralel omp sections directive

  finish=clock();

  cout<<"Nos tardamos "<<finish-start<< " Ciclos de cpu "<<endl; 
  
  return 0;
 
}

