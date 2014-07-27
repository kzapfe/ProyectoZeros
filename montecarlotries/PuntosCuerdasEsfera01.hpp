//Veamos que podemos hacer con esto
//Tratemos de encontrar los Blind Spots Subplanckianos
//con esta chiva.


//No necesitas incluir nada que ya se haya incluido antes
#include "RutinasExternas02.hpp"

#include "intervalitos01.hpp"
#include "ComparaWeylenradio01.hpp"

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
  

  //***************************************************************************
  // Constantes
  //Solo piensa en superposicion de ondas medio aleatoreas
  const int resolucion=10; //La resolucion es una propiedad de la rutina.
  const double epsilon=0.01;
  const double tolzero=0.000000;




  
  //***************************************************************************
  // Esta es la parte en donde ponemos los nombres de los archivos y los abrimos
  //***************************************************************************

  //Deberias poner esto en una rutina externa
  //Es una patada en los huevos

  std::ostringstream escupechord;
  std::ostringstream escupecondini;
  std::ostringstream escuperadial;
  std::ostringstream escupezeros;

  escupechord<<muestreo<<"_ChordsSphereT.dat"<<std::ends;
  escupecondini<<muestreo<<"_condiniSphere.dat"<<std::ends;
  escuperadial<<muestreo<<"_ChordsSphereHomogenT.dat"<<std::ends;
  escupezeros<<muestreo<<"_ZerosSphereTzolkin.dat"<<std::ends;

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


  //***************************************************************************
  // Poblamos la esfera
  //***************************************************************************
  //Bueno, si el sistema es caotico 
  //tiene al menos dos grados de libertad
  simplectic * x;
  simplectic * y;
0  

  
  x=new simplectic[muestreo];
  y=new simplectic[muestreo];
  mat bolallena;
  
  bolallena=Populate4BallShell(nivel, muestreo);
  

  for(int i=0;i<muestreo;i++){
    x[i].q=bolallena(i,0);
    x[i].p=bolallena(i,1);
    y[i].q=bolallena(i,2);
    y[i].p=bolallena(i,3);
  };
  
  cout<<" Acabamos de popular la capa de energia"
      << " con "<<muestreo<<" puntos"<<endl;
  
  int check;
  check=min(muestreo,2756);

  
  for(int i=0;i<check;i++){
    //For check, we also print the magnitude
    condinispace<<x[i].q<<"\t"<<x[i].p<<"\t"
		<<y[i].q<<"\t"<<y[i].p<<"\t"
		<<x[i].q*x[i].q+x[i].p*x[i].p+
      y[i].q*y[i].q+y[i].p*x[i].p
		<<endl;
  };

  
  //***************************************************************************
  
  // Variables auxiliares para encontrar el zero.
  double bigrho=1.0, smallrho=0.00000;
  double theta1, theta2, theta3,rho; 
  intervalito InterVal;

  for(int i=0;i<resolucion;i++){
    cout<<" pasando al theta1 \t" << i <<endl;
    
    for(int k=0;k<resolucion;k++){	  
      for(int l=-resolucion;l<resolucion;l++){
	gsl_complex Weyl;
	double weylreal=0.00;
	double weylimag=0.00;
	simplectic xhi(0.0,0.0), mu(0.000,0.000);

	InterVal.menor=smallrho;
	InterVal.mayor=bigrho;
	
	
	  //estos corren de zero a pi
	  theta1=pi*(double)i/(double)resolucion;
	  theta2=pi*(double)k/(double)resolucion;
	  
	  //este corre de -pi a pi
	  theta3=pi*(double)l/(double)resolucion;
	       
	  //------------------------------------------------
	  //Evaluamos la funcion en los dos extremos
	  rho=InterVal.menor;
  
	  mu.q= cos(theta1)*rho;
	  mu.p= sin(theta1)*cos(theta2)*rho;
	  xhi.q=sin(theta1)*sin(theta2)*cos(theta3)*rho;
	  xhi.p=sin(theta1)*sin(theta2)*sin(theta3)*rho; 
	  
	  Weyl=WeylFunction2D(muestreo, x, y, mu, xhi);
	  //weylreal=GSL_REAL(Weyl);
	  //weylimag=GSL_IMAG(Weyl);

	  chordradialspace<<rho<<"\t"<<theta1<<"\t"
	  		  <<theta2<<"\t"<<theta3<<"\t"
	  		  <<weylreal<<"\t"<<weylimag<<endl;
	  chordspace<<mu.q<<"\t"<<mu.p<<"\t"
		    <<xhi.q<<"\t"<<xhi.p<<"\t"
		    <<weylreal<<"\t"<<weylimag<<endl;


	  rho=InterVal.mayor;
  
	  cout<<rho<< " THIS IS DAMN RHO "<<endl;
	  mu.q= cos(theta1)*rho;
	  mu.p= sin(theta1)*cos(theta2)*rho;
	  xhi.q=sin(theta1)*sin(theta2)*cos(theta3)*rho;
	  xhi.p=sin(theta1)*sin(theta2)*sin(theta3)*rho; 
	  
	  Weyl=WeylFunction2D(muestreo, x, y, mu, xhi);
	  weylreal=GSL_REAL(Weyl);
	  weylimag=GSL_IMAG(Weyl);
	  
	  chordradialspace<<rho<<"\t"<<theta1<<"\t"
			  <<theta2<<"\t"<<theta3<<"\t"
			  <<weylreal<<"\t"<<weylimag<<endl;
	  
	  chordspace<<mu.q<<"\t"<<mu.p<<"\t"
		    <<xhi.q<<"\t"<<xhi.p<<"\t"
		    <<weylreal<<"\t"<<weylimag<<endl;
	  //-----------------------------------------------

	  
	  InterVal.cambia=1;
	 
	  if(InterVal.cambia)cout<<InterVal.cambia<<"\t"<<i<<endl;
	  
	       
	  if(InterVal.cambia){
	    cout<<"cambio!"<<endl;
	    //Aqui buscamos los Zeros
	    //Tenemos que diezmar los datos
	    // Si no el puto archivo diezmalineas 
	    // queda monstruosamente grande
	    const int profundidad=5;
	    int divisiones=5;
	    int iteracion=0;
	    intervalito auxinterval=InterVal;
	    double longitud=InterVal.mayor-InterVal.menor;
	    
	    while((iteracion<profundidad)&&(longitud>epsilon)){
	      
	      intervalito *subintervalos;
	      
	      subintervalos=Particion(divisiones, auxinterval);
	      cout<<iteracion<< " iteracion" <<endl;

	      for(int pieza=0; pieza<divisiones; pieza++){
		
		subintervalos[pieza].cambia=ComparaWeylRadial2D
		  (subintervalos[pieza], theta1,  theta2,  theta3, 
		   x, y, muestreo);
		if(subintervalos[pieza].cambia){
		  auxinterval=subintervalos[pieza];
		  iteracion++;
		  pieza=divisiones+1;
		};
	      };
	    
	    };
	    
	    zeroline<<auxinterval.menor<<"\t"<<auxinterval.mayor<<"\t"
		    <<theta1<<"\t"<<theta2<<"\t"<<theta3<<"\t"
		    <<endl;
      
	    

	    
	  };
	  

	  
      };
      //chordspace<<endl;  
      //chordradialspace<<endl;  	
    };
  zeroline<<endl;
  zeroline<<endl;
  };
  


  
  


  chordspace.close();
  condinispace.close();
  chordradialspace.close();
  zeroline.close();

  return;


};
