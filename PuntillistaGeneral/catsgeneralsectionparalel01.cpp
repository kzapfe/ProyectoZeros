//Probemos hacer gato estados cuadridimensionales.
//General version, carga los datos de los centros a partir de un archivo externo


#include <iostream>
#include <fstream>
#include <vector>
#include <armadillo>

#include <ctime> //para ver cuanto se tarda en cada punto
#include <omp.h> //para ver cuanto se tarda en cada punto

//Ultra Global Variables

const int resol=100;
double magnituddeinteres=2.7;

#include "PopulateSphericalShells01.hpp"
#include "CortesGatosBidimensionales01.hpp"
#include "CatStates01.hpp"

using namespace std;
using namespace arma;

/********************************************************************/

int main(int argc, char* argv[]){
  //No mames que bien trabajas despues de hacer ejercicio 
  // Y SIN INTERNET. Recuerda eso SIEMPRE.
  int maximumgauss=1;
  vector   <simplectic> x,y;
  string NomdeCentros;

  //El numero de estados coherentes. El numero de gatos es 
  // la combinacion de pares posibles.

  //Los centros, otra vez para variar.
  mat centros;
  centros.load(NomdeCentros);
  maximumgauss=centros.n_rows;

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
      testeando="CorteExactGaussQchico.dat";

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
	   testeando="CorteExactGaussXchico.dat";

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
	   testeando="CorteExactGaussYchico.dat";

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
	   testeando="CorteExactGaussPchico.dat";

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
      testeando="CorteExactGaussCuerdasQchico.dat";

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
	   testeando="CorteExactGaussCuerdasMuchico.dat";

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
	   testeando="CorteExactGaussCuerdasXhichico.dat";

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
	   testeando="CorteExactGaussCuerdasPchico.dat";

	   WeylSection(CentroX, CentroY, GatosX, GatosY, 
		    maximumgauss, cuentainterferencias,
		    qx, qy, testeando);

   
      
    } //Cierra el ultimo  pragma omp section de cuerdas

  } //Cierra el pragma paralel omp sections directive

  finish=clock();

  cout<<"Nos tardamos "<<finish-start<< " Ciclos de cpu "<<endl; 
  
  return 0;
 
}

