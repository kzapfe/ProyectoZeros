//Veamos que podemos hacer con esto
//Tratemos de encontrar los Blind Spots Subplanckianos
//con esta chiva.


//No necesitas incluir nada que ya se haya incluido antes
//Mais uma vez, esta vez con capacidad de encontar todos los
//zeros en el intervalo iniial.


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
  const int resolucion1=25; //La resolucion es una propiedad de la rutina.
  const int resolucion2=5; //La resolucion es una propiedad de la rutina.
  const int resolucion3=1; //La resolucion es una propiedad de la rutina. 
  //***************************************************************************
  // Esta es la parte en donde ponemos los nombres de los archivos y los abrimos
  //***************************************************************************

  //Deberias poner esto en una rutina externa
  //Es una patada en los huevos

  std::ostringstream escupechord;
  std::ostringstream escupecondini;
  std::ostringstream escuperadial;
  std::ostringstream escupezeros;

  escupechord<<muestreo<<"_VerificarFuncion.dat"<<std::ends;
  escupecondini<<muestreo<<"_condiniSphere.dat"<<std::ends;
  escuperadial<<muestreo<<"_ChordsSphereHomogenT.dat"<<std::ends;
  escupezeros<<muestreo<<"_ZerosVarios3Tz.dat"<<std::ends;

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
  for(int m=0; m<8;m++){ 
    // Variables auxiliares para encontrar el zero.

    double bigrho=0.01+0.15*(double)(m+1), smallrho=0.01+0.15*(double)m;
    double theta1, theta2, theta3,rho; 
    intervalito InterVal;

    InterVal.menor=smallrho;
    InterVal.mayor=bigrho;
    cout<<m<<"th Interval "<< smallrho << "\t" <<bigrho<<endl;
    
    for(int i=0;i<resolucion1;i++){
      //for(int i=5;i<6;i++){
      
      cout<<" pasando al theta1 \t" << i <<endl;
    
      for(int k=0;k<resolucion2;k++){	  
	for(int l=-resolucion3;l<resolucion3;l++){
	  
	  //Temporary Thest
	  //for(int k=5;k<6;k++){	  
	  // for(int l=0;l<1;l++){
	  
	  
	  gsl_complex Weyl;
	  double weylreal=0.00;
	  double weylimag=0.00;
	  simplectic xhi(0.0,0.0), mu(0.000,0.000);
	  
	  
	  //estos corren de zero a pi
	  theta1=pi*(double)i/(double)resolucion1;
	  theta2=pi*(double)k/(double)resolucion2;
	  
	  //este corre de -pi a pi
	  theta3=pi*(double)l/(double)resolucion3;
	  
	  //Aqui buscamos los Zeros
	  //Tenemos que diezmar los datos
	  // Si no el puto archivo diezmalineas 
	  // queda monstruosamente grande
	  const int profundidad=6;
	  int divisiones=7;
	  int iteracion=0;
	  intervalito auxinterval=InterVal;
	  //double longitud=InterVal.mayor-InterVal.menor;
	  int pieza;
	  
	  intervalito *subintervalos;
	  subintervalos=Particion(divisiones, auxinterval);
	  

	  while(iteracion<profundidad){
	    
	    /*
	      cout<<"Este intervalo va de "<<auxinterval.menor
	      << " hasta "<<auxinterval.mayor
	      <<", vamos a dividirlo en piezas"<<endl;
	    */ 
	    
	    
	    subintervalos=Particion(divisiones, auxinterval);
	    pieza=0;
	    
	    
	    
	    while(pieza<divisiones){
	    /*
	      cout<<"Este subintervalo va de "<<subintervalos[pieza].menor
	      << " hasta "<<subintervalos[pieza].mayor<<endl;
	    */
	      //Temporarly this merde here, cause if not, doesnt work.
	      
	      bool result;
	      double rho;
	      simplectic xhi, mu;
	      int bigrhoWeylSign=0;
	      int smallrhoWeylSign=0;
	      
	      double weylreal=0.0;
	      double weylimag=0.0;
	      
	      
	     //**************************************************
	     //Get the diabolical Weyl Function
	     //for this mu and xhi
	     rho=subintervalos[pieza].menor;
	     
	     mu.q= cos(theta1)*rho;
	     mu.p= sin(theta1)*cos(theta2)*rho;
	     xhi.q=sin(theta1)*sin(theta2)*cos(theta3)*rho;
	     xhi.p=sin(theta1)*sin(theta2)*sin(theta3)*rho;
  
	     //Temporarly this merde here, cause if not, doesnt work.
	     
	     for(int k=0; k<muestreo; k++){
	       //Aqui hacemos la transformada de Fourier
	       // Que nos da la Funcion de Weyl
	       
	       weylreal+=cos((x[k].simplecticproduct(mu)+
			      y[k].simplecticproduct(xhi))/hbarra);
	       
	       weylimag+=sin((x[k].simplecticproduct(mu)+
			      y[k].simplecticproduct(xhi))/hbarra);
    
	     };
	     
	     
	     
	     weylreal=weylreal/(hbarra*pi*muestreo);
	     weylimag=weylimag/(hbarra*pi*muestreo);
	     
	     //Weyl=WeylFunction2D(muestreo, x, y, mu, xhi);
	     
	     //weylreal=GSL_REAL(Weyl);
	     //weylimag=GSL_IMAG(Weyl);
	     

	     
	     smallrhoWeylSign=signum(weylreal);
	     chordspace<<rho<<"\t"<<weylreal<<
	       "\t"<<iteracion<<endl; 
	      

	     //**************************************************
	     
	     //**************************************************
	     //Get the diabolical Weyl Function
	     //for this mu and xhi
	     rho=subintervalos[pieza].mayor;
	     
	     mu.q= cos(theta1)*rho;
	     mu.p= sin(theta1)*cos(theta2)*rho;
	     xhi.q=sin(theta1)*sin(theta2)*cos(theta3)*rho;
	     xhi.p=sin(theta1)*sin(theta2)*sin(theta3)*rho;
	     
	     for(int k=0; k<muestreo; k++){
	       //Aqui hacemos la transformada de Fourier
	       // Que nos da la Funcion de Weyl
	       
	       weylreal+=cos((x[k].simplecticproduct(mu)+
		   y[k].simplecticproduct(xhi))/hbarra);
	       
	       weylimag+=sin((x[k].simplecticproduct(mu)+
		   y[k].simplecticproduct(xhi))/hbarra);
	       
	     };
  
  
  
	     weylreal=weylreal/(hbarra*pi*muestreo);
	     weylimag=weylimag/(hbarra*pi*muestreo);
	  
	     chordspace<<rho<<"\t"<<weylreal<<
		       "\t"<<iteracion<<endl; 
	     chordspace<<endl;

	     bigrhoWeylSign=signum(weylreal);
  
	     //**************************************************	  
	     
	     result=(bigrhoWeylSign!=smallrhoWeylSign);
	     subintervalos[pieza].cambia=result;

	     //FIN DA MERDE

	     /*
	    subintervalos[pieza].cambia=ComparaWeylRadial2D
	      (subintervalos[pieza], theta1,  theta2,  theta3, 
	       x, y, muestreo);
	     */



	    
	      if(subintervalos[pieza].cambia){
		auxinterval=subintervalos[pieza];
		auxinterval.cambia=true;
		
		/*cout<< "Hubo cambio aqui,pieza "<<pieza<<
		  " en la iteracion: "
		    <<iteracion<<endl;
		*/

		pieza=divisiones+1;
		iteracion++;
		
	      }else{
		pieza++;
		//cout<< "No cambio aqui, check next pieza: "<<pieza<<endl;
		if(pieza==divisiones){
		  /*cout<<"No hubo cambios algo esta mal en la iteracion "
		      <<iteracion<<endl;
		      cout<<endl;*/
		  iteracion++;
		};
	      }; 
	      
	  };
	    
	};
	
	zeroline<<auxinterval.menor<<"\t"<<theta1<<"\t"
		<<theta2<<"\t"
		<<theta3<<"\t"
		<<iteracion<<"\t"
		<<auxinterval.cambia
		<<endl;
	    

	
	
	
	
	
	  

	  
      };
      //chordspace<<endl;  
      //chordradialspace<<endl;  	
      };
    zeroline<<endl;
    zeroline<<endl;
  };
  
  };

  
  


  chordspace.close();
  condinispace.close();
  chordradialspace.close();
  zeroline.close();

  return;


};
