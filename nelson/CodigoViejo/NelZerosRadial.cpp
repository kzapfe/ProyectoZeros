/*
Este programa dibuja la funcion de Cuerdas para 
planos mu=cte, en subdominio xhi. 
Try to make things cleaner.
 */

/*esta modificacion va a buscar los zeros radialmente
 */

#include "RutinasNelson02.hpp"
#include "intervalito02.hpp"
//Los parametros globales en su propio archivo.
#include "ParametrosGlobales.hpp"

using namespace std;



int main(){  
  //inicializar el semillador
  srand(12389);  
  
  //Abrimos los archivos necesarios_-----------------------------
  
  std::string hazstring;
  
  std::ostringstream escupezreal;
  escupezreal<<muestreo<<"_"<<nivel<<
    "_ZeroNelson-Real-Rad-0-0.dat"<<std::ends;
  hazstring=escupezreal.str();
  const char *nombrezeroreal=hazstring.c_str();
  ofstream ZeroReal;
  ZeroReal.open(nombrezeroreal);
  
  std::ostringstream escupezimag;
  escupezimag<<muestreo<<"_"<<nivel<<
    "_ZeroNelson-Imag-Rad-0-0.dat"<<std::ends;
  hazstring=escupezimag.str();
  const char *nombrezeroimag=hazstring.c_str();
  ofstream ZeroImag;
  ZeroImag.open(nombrezeroimag);
  
  std::ostringstream escupemedias;
  escupemedias<<muestreo<<"_"<<nivel<<
    "_MediasWignerZero0-0.dat"<<std::ends;
  hazstring=escupemedias.str();
  const char *nombremedwig=hazstring.c_str();
  ofstream MediasWigner;
  MediasWigner.open(nombremedwig);
  
  
  ofstream testa;
  testa.open("Datosparachecar.dat");
  
  
  //----------------------------------------------------
  
  //Poblar la camada de energia
  mat DiracDeltas;
  DiracDeltas=PopulateNelson(nivel, muestreo);
    
  cout<<" Acabamos de popular la capa de energia"
      << " con "<<muestreo<<" puntos"<<endl;
  
  mat meanvalues;
  
  meanvalues=mean(DiracDeltas);
  
  MediasWigner<<"Primer momento x_q: "<< meanvalues(0)<<endl;
  MediasWigner<<"Primer momento x_p: "<< meanvalues(1)<<endl;
  MediasWigner<<"Primer momento y_q: "<< meanvalues(2)<<endl;
  MediasWigner<<"Primer momento y_p: "<< meanvalues(3)<<endl;
  MediasWigner.close();


  simplectic * x;
  simplectic * y;
  
  x=new simplectic[muestreo];
  y=new simplectic[muestreo];
  
  
  for(int i=0;i<muestreo;i++){
    x[i].q=DiracDeltas(i,0);
    x[i].p=DiracDeltas(i,1);
    y[i].q=DiracDeltas(i,2);
    y[i].p=DiracDeltas(i,3);
  };
  
  

  //----------------------------------------------------
  
  //********HERE BEGINS THE PARTY !!!******************
  //----------------------------------------------------

  intervalito LineaRadial;
  LineaRadial.menor=0;
  LineaRadial.mayor=extrem;
  
  intervalito *subintervalo;
  subintervalo=Particion(divisiones, LineaRadial);
  
  double cuerdasrealchico, cuerdasrealgrande;
  simplectic xhi(0.0,0.0);
  
  
  //--------------------------------------------
  //primero, localizamos los subintervalos
  //con cambio de signo
  //no mas de divisiones, por favor

  intervalito storeinterval[divisiones];
  int counter=0;
  int const semires=20;

  for(int j=0;j<resolucion; j++){
    double theta=pi*j/(2.0*(double)resolucion);
    for(int i=0; i<divisiones; i++){    
      double rho=0.00;
	 
    
      rho=subintervalo[i].menor;
      xhi.q=rho*cos(theta);
      xhi.p=rho*sin(theta);
      cuerdasrealchico=WeylReal(muestreo, mu, xhi, x, y);
      
      
      rho=subintervalo[i].mayor;
      xhi.q=rho*cos(theta);
      xhi.p=rho*sin(theta);
      cuerdasrealgrande=WeylReal(muestreo, mu, xhi, x, y);
      
      subintervalo[i].cambia=
	(signum(cuerdasrealgrande)!=signum(cuerdasrealchico));
      
      if(subintervalo[i].cambia){
	storeinterval[counter]=subintervalo[i];
	counter++;
      };
      
      cout<<"Hasta ahora no ha pasado nada "<<theta<<endl;
    };

    //--------------------------------------------------


    //Now to get deeper on each such interval
  
  
    for(int i=0; i<counter; i++){ 
      intervalito auxinterval=storeinterval[i];
      bool checando=false;
      
      
    //Now the real search----6660----
      int deep=0;
      int pieza=0;
      double valorweyl;
      
      subintervalo=Particion(divisiones, auxinterval);
      
      while((pieza<divisiones)&&(deep<maxprof)){

	
	double rho=0.00;
	    
	rho=subintervalo[i].menor;
	xhi.q=rho*cos(theta);
	xhi.p=rho*sin(theta);
	cuerdasrealchico=WeylReal(muestreo, mu, xhi, x, y);
	
	
	rho=subintervalo[i].mayor;
	xhi.q=rho*cos(theta);
	xhi.p=rho*sin(theta);
	cuerdasrealgrande=WeylReal(muestreo, mu, xhi, x, y);
      

	subintervalo[pieza].cambia=
	  (signum(cuerdasrealgrande)!=signum(cuerdasrealchico));
	if(subintervalo[pieza].cambia){	  
	  if(deep<maxprof){
	    auxinterval=subintervalo[pieza];
	    subintervalo=Particion(divisiones, auxinterval);
	    deep++;	    
	    xhi.p=auxinterval.menor;
	    valorweyl=cuerdasrealchico;
	    checando=true;
	    pieza=0;
	  }else{pieza=divisiones+1;};	  
	}else{
	  pieza++;	  
	} ;
	
      };
      
      if(checando)ZeroReal<<xhi.q<<"\t"<<xhi.p<<"\t"<<valorweyl<<endl;
    
  
      
    };//Este termina el counter sobre los posibles zeros.
  };//este cierra sobre las barras verticales


  return 0;

    
}
