/*
Este programa dibuja la funcion de Cuerdas para 
planos mu=cte, en subdominio xhi. 


 */
#include <armadillo>
#include "simplectic01.hpp"
#include "ParametrosGlobales.hpp"
#include "RutinasNelson02.hpp"


using namespace std;



int main(){
  
  //inicializar el semillador
  srand(12389); 
  cout<<"la seccion tiene q y p: "<<endl;
  cout<<mu.q<<"\t"<<mu.p<<endl;
  //Abrimos los archivos necesarios

  std::ostringstream escupeweyl;
  std::string hazstring;

  escupeweyl<<"TresScars_294"<<"_"<<"08138"<<
    "_Weyl0-0.dat"<<std::ends;
  hazstring=escupeweyl.str();
  const char *nombreweyl=hazstring.c_str();
  ofstream WeylSeccion;
  WeylSeccion.open(nombreweyl);

 

  std::ostringstream escupemedias;

  escupemedias<<"TresScars_294"<<"_"<<"08138"<<
    "_Medias0-0.dat"<<std::ends;
  hazstring=escupemedias.str();
  const char *nombremedwig=hazstring.c_str();
  ofstream MediasWigner;
  MediasWigner.open(nombremedwig);


  ofstream Centros;
  Centros.open("CentrosTresScars.dat");



  //----------------------------------------------------
  //Poblar la camada de energia
  mat DiracDeltas;
  mat auxiliar1,auxiliar2,auxiliar3;

  double Proporcion;
  double Tau1,Tau2,Tau3;

  //Periodo de la orbita vertical
  Tau1=2.0*pi;

  //Periodo de la orbita que comienza en 0.48854
  Tau2=11.6600;
  
  //Periodo de la orbita que comienza en 0.723954
  Tau3=6.4867;


  
  Proporcion=(double)muestreo/(1/Tau1+1/Tau2+1/Tau3);
  
  int muestreo1, muestreo2, muestreo3;

  muestreo1=Proporcion/Tau1;
  muestreo2=Proporcion/Tau2;
  muestreo3=muestreo-muestreo1-muestreo2;

  cout<<"Este es muestreo2: "<<muestreo2<<endl;

  auxiliar1=PopulateCircle(nivel, muestreo1);  
  auxiliar2= PopulatefromDrive(muestreo2, "TrayectoriaPeriodica02.dat");
  auxiliar3= PopulatefromDrive(muestreo3, "TrayectoriaPeriodica03.dat");
  

  DiracDeltas=zeros<mat>(muestreo,4);


  DiracDeltas.submat(0,0,muestreo1-1,3)=auxiliar1;
  DiracDeltas.submat(muestreo1,0,muestreo1+muestreo2-1,3)=auxiliar2;
  DiracDeltas.submat(muestreo1+muestreo2,0,muestreo-1,3)=auxiliar3;

  
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
  

  for(int i=0; i<muestreo;i++){
    //Cheq that everything is as it should
    Centros<<x[i].q<<"\t"<<x[i].p<<"\t"
	   <<y[i].q<<"\t"<<y[i].p<<"\t"
		<<NelsonEnergy(DiracDeltas.row(i))<<endl;
  };


  cout<<" Acabamos de poblar la cicatriz mas pequenha"
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

	  

	  weylreal+=cos(-(x[k].simplecticproduct(mu)+
			 y[k].simplecticproduct(xhi))/hbarra);
	  
	  weylimag+=sin(-(x[k].simplecticproduct(mu)+
			 y[k].simplecticproduct(xhi))/hbarra);
	  
	};
	
	     
       	weylreal=weylreal/(hbarra*pi*muestreo);
	weylimag=weylimag/(hbarra*pi*muestreo);
	
	WeylSeccion<<xhi.q<<"\t"<<xhi.p<<
	  "\t"<<weylreal<<"\t"<<weylimag<<endl;
	
		
	
    };
    
    
    WeylSeccion<<endl;
    
    
    cout<<"Vamos en la linea i="<<i<<endl;
  
  };

  
  WeylSeccion.close();
  
 
  return 0;
  
}
