/*
Este programa dibuja la funcion de Cuerdas para 
planos mu=cte, en subdominio xhi. 
 */

#include "RutinasNelson01.hpp"
#include "intervalito02.hpp"

using namespace std;

//Esto es pleonasmo, pero para que te acuerdes.

int main(){
  
  //inicializar el semillador
  srand(12389);

  //Parametros Globales
  int muestreo=500000;
  int nivel=10;
  int resolucion=100;
  double extrem=1.500;
  double tol=0.000500;

  //seciones mu constante!
  simplectic  mu(0.000,0.000);
  mu.q= 0.00;
  mu.p= 0.00;
  
  //Abrimos los archivos necesarios

  std::ostringstream escupeweyl;
  std::string hazstring;

  escupeweyl<<muestreo<<"_"<<nivel<<
    "_NelsonW0-0.dat"<<std::ends;
  hazstring=escupeweyl.str();
  const char *nombreweyl=hazstring.c_str();
  ofstream WeylSeccion;
  WeylSeccion.open(nombreweyl);



  std::ostringstream escupezreal;

  escupezreal<<muestreo<<"_"<<nivel<<
    "_ZeroRealNelson0-0.dat"<<std::ends;
  hazstring=escupezreal.str();
  const char *nombrezeroreal=hazstring.c_str();
  ofstream ZeroReal;
  ZeroReal.open(nombrezeroreal);


  std::ostringstream escupezimag;

  escupezimag<<muestreo<<"_"<<nivel<<
    "_ZeroImagNelson0-0.dat"<<std::ends;
  hazstring=escupezimag.str();
  const char *nombrezeroimag=hazstring.c_str();
  ofstream ZeroImag;
  ZeroImag.open(nombrezeroimag);






  std::ostringstream escupemedias;

  escupemedias<<muestreo<<"_"<<nivel<<
    "_MediasWigner0-0.dat"<<std::ends;
  hazstring=escupemedias.str();
  const char *nombremedwig=hazstring.c_str();
  ofstream MediasWigner;
  MediasWigner.open(nombremedwig);




  //----------------------------------------------------
  //Poblar la camada de energia
  mat DiracDeltas;
  DiracDeltas=PopulateNelson(nivel, muestreo);
  DiracDeltas.save("DiracNelson.dat", raw_ascii);
  
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
  
  cout<<" Acabamos de popular la capa de energia"
      << " con "<<muestreo<<" puntos"<<endl;
  
  mat meanvalues;
  
  meanvalues=mean(DiracDeltas);
  
  MediasWigner<<"Primer momento x_q: "<< meanvalues(0)<<endl;
  MediasWigner<<"Primer momento x_p: "<< meanvalues(1)<<endl;
  MediasWigner<<"Primer momento y_q: "<< meanvalues(2)<<endl;
  MediasWigner<<"Primer momento y_p: "<< meanvalues(3)<<endl;
  MediasWigner.close();


  


  //----------------------------------------------------
  
  //********HERE BEGINS THE PARTY !!!******************
  //----------------------------------------------------


      
    

  for(int i=-resolucion;i<resolucion;i++){
    for(int k=-resolucion;k<resolucion;k++){
	   
	 
	double weylreal=0.00;
	double weylimag=0.00;
	simplectic xhi(0.0,0.0);
	  
	     
	
	xhi.q=(double)i*extrem/(double)resolucion;
	xhi.p=(double)k*extrem/(double)resolucion;
	
	
	for(int k=0; k<muestreo; k++){
	  //Aqui hacemos la transformada de Fourier
	  // Que nos da la Funcion de Weyl
	  /// guuuey, es xhi^x, no al revez....
	  // La parte real es simetrica, pero la otra...
	  
	  weylreal+=cos(-(x[k].simplecticproduct(mu)+
			 y[k].simplecticproduct(xhi))/hbarra);
	  
	  weylimag+=sin(-(x[k].simplecticproduct(mu)+
			 y[k].simplecticproduct(xhi))/hbarra);
	  
	};
	
	     
       	weylreal=weylreal/(hbarra*pi*muestreo);
	weylimag=weylimag/(hbarra*pi*muestreo);
	
	WeylSeccion<<xhi.q<<"\t"<<xhi.p<<
	  "\t"<<weylreal<<"\t"<<weylimag<<endl;
	
	if(abs(weylreal)<tol){
	  ZeroReal<<xhi.q<<"\t"<<xhi.p<<
	  "\t"<<weylreal<<"\t"<<weylimag<<endl;
	};

	
	if(abs(weylimag)<tol){
	  ZeroImag<<xhi.q<<"\t"<<xhi.p<<
	  "\t"<<weylreal<<"\t"<<weylimag<<endl;
	};

		
	
    };

    WeylSeccion<<endl;
    ZeroReal<<endl;
    ZeroImag<<endl;
    
    cout<<"Vamos en la linea i="<<i<<endl;
  
  };

  
  WeylSeccion.close();
  ZeroReal.close();
  ZeroImag.close();
 
  return 0;
  
}
