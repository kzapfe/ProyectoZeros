/*
Este programa dibuja la funcion de Cuerdas para 
planos mu=cte, en subdominio xhi. 
 */

#include "RutinasNelson02.hpp"
#include "intervalito02.hpp"

using namespace std;

//Esto es pleonasmo, pero para que te acuerdes.

int main(){
  
  //inicializar el semillador
  srand(12389);

  //Parametros Globales
  int muestreo=5000;
  int nivel=10;
  int resolucion=250;
  double extrem=1.500;
  

  //seciones mu constante!
  simplectic  mu(0.000,0.000);
  mu.q= 0.00;
  mu.p= 0.00;
  
  //Abrimos los archivos necesarios


  std::string hazstring;


  std::ostringstream escupezreal;

  escupezreal<<muestreo<<"_"<<nivel<<
    "_ZeroNelson-Real-V-0-0.dat"<<std::ends;
  hazstring=escupezreal.str();
  const char *nombrezeroreal=hazstring.c_str();
  ofstream ZeroReal;
  ZeroReal.open(nombrezeroreal);


  std::ostringstream escupezimag;

  escupezimag<<muestreo<<"_"<<nivel<<
    "_ZeroNelson-Imag-V-0-0.dat"<<std::ends;
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
  testa.open("Datosaloguey.dat");
    
    
  //----------------------------------------------------
  //Poblar la camada de energia
    mat DiracDeltas;
  DiracDeltas=PopulateNelson(nivel, muestreo);
  
  
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

  const int divisiones=7;
  const int maxprof=7;
   

  //REAL PARTY
  for(int i=-resolucion; i<resolucion; i++){
   
    //estas se reinicialican en cada columna
    intervalito InterVal;
    intervalito auxinterval;
    
    
    simplectic xhi(0.0,0.0);
    double weylreal;


    InterVal.menor=-extrem;
    InterVal.mayor=extrem;
    auxinterval=InterVal;
    

    intervalito *subintervalo;


    int pieza, iteracion=0;


    while(iteracion<maxprof){

      pieza=0;
      subintervalo=Particion(divisiones, auxinterval);
     

      while(pieza<divisiones){  
	
       
	double weylrealmayor, weylrealmenor;

       
      
	//rastreo normal en la direccion horizontal
	xhi.q=(double)i*extrem/(double)resolucion;
	
	//rastreo por signo en la vertical
	
       
	
	//Pieza Menor__________________________________________
	xhi.p=subintervalo[pieza].menor;	  
	
	weylreal=0.00;
		
	for(int k=0; k<muestreo; k++){
	  //Evaluando extremo inferior del intervalo
	  
	  weylreal+=cos(-(x[k].simplecticproduct(mu)+
			  y[k].simplecticproduct(xhi))/hbarra);
	  
	};

	weylrealmenor=weylreal/(hbarra*pi*muestreo);

	testa<<xhi.q<<"\t"<<xhi.p<<"\t"<<weylrealmenor<<
	  "\t"<<pieza<<"\t"<<iteracion<<endl;

	//Pieza Mayor_______________________________________
	xhi.p=subintervalo[pieza].mayor;
	weylreal=0.00;

	for(int k=0; k<muestreo; k++){
	  //Evaluando extremo superior del intervalo
	  weylreal+=cos(-(x[k].simplecticproduct(mu)+
			  y[k].simplecticproduct(xhi))/hbarra);	  
	};
	weylrealmayor=weylreal/(hbarra*pi*muestreo);


	//__________________________________________________

	
	testa<<xhi.q<<"\t"<<xhi.p<<"\t"<<weylrealmayor<<
	  "\t"<<pieza<<"\t"<<iteracion<<endl;
	testa<<endl;

	subintervalo[pieza].cambia=
	    (signum(weylrealmayor)!=signum(weylrealmenor));
	
	  
	  if(subintervalo[pieza].cambia){
	    auxinterval=subintervalo[pieza];
	    auxinterval.cambia=true; 
	    
	    if(iteracion==maxprof-1){
	      ZeroReal<<xhi.q<<"\t"<<xhi.p<<
		"\t"<<weylreal<<"\t"<<endl;};

	    iteracion++;
	    pieza=divisiones;
	  }else if(pieza>=divisiones-1){		
	    //   cout<<"Next ray, pleze"<<endl;
	    iteracion=maxprof+1;   
	  }else{pieza++;};
	  
      };    
    };


    cout<<"Vamos en la linea i="<<i<<endl;
	
  };

  ZeroReal.close();

  //ENDS REAL PARTY

  
  //Imaginary PARTY
  for(int i=-resolucion; i<resolucion; i++){
   
    //estas se reinicialican en cada columna
    intervalito InterVal;
    intervalito auxinterval;
    
    
    simplectic xhi(0.0,0.0);
    double weylimag;


    InterVal.menor=-extrem;
    InterVal.mayor=extrem;
    auxinterval=InterVal;
    

    intervalito *subintervalo;


    int pieza, iteracion=0;


    while(iteracion<maxprof){

      pieza=0;
      subintervalo=Particion(divisiones, auxinterval);
      cout<<"la la la estamos locos"<<endl;

      while(pieza<divisiones){  
	
       
	double weylimagmayor, weylimagmenor;
      
	//rastreo normal en la direccion horizontal
	xhi.q=(double)i*extrem/(double)resolucion;
	
	//rastreo por signo en la vertical
	
       
  	//Pieza Menor__________________________________________
	xhi.p=subintervalo[pieza].menor;	  
	weylimag=0.00;		
	for(int k=0; k<muestreo; k++){
	  //Evaluando extremo inferior del intervalo
	  weylimag+=sin(-(x[k].simplecticproduct(mu)+
			  y[k].simplecticproduct(xhi))/hbarra);
	  
	};

	weylimagmenor=weylimag/(hbarra*pi*muestreo);

	testa<<xhi.q<<"\t"<<xhi.p<<"\t"<<weylimagmenor<<
	  "\t"<<pieza<<"\t"<<iteracion<<endl;

	//Pieza Mayor_______________________________________
	xhi.p=subintervalo[pieza].mayor;
	weylimag=0.00;

	for(int k=0; k<muestreo; k++){
	  //Evaluando extremo superior del intervalo
	  weylimag+=sin(-(x[k].simplecticproduct(mu)+
			  y[k].simplecticproduct(xhi))/hbarra);	  
	};
	weylimagmayor=weylimag/(hbarra*pi*muestreo);

	//__________________________________________________

	
	testa<<xhi.q<<"\t"<<xhi.p<<"\t"<<weylimagmayor<<
	  "\t"<<pieza<<"\t"<<iteracion<<endl;
	testa<<endl;

	subintervalo[pieza].cambia=
	    (signum(weylimagmayor)!=signum(weylimagmenor));

	
	  
	  if(subintervalo[pieza].cambia){
	    auxinterval=subintervalo[pieza];
	    auxinterval.cambia=true; 
	    
	    if(iteracion==maxprof-1){
	      ZeroImag<<xhi.q<<"\t"<<xhi.p<<
		"\t"<<"\t"<<weylimag<<endl;
	    };

	    iteracion++;
	    pieza=divisiones;
	  }else if(pieza>=divisiones-1){		
	    //   cout<<"Next ray, pleze"<<endl;
	    iteracion=maxprof+1;   
	  }else{pieza++;};
	  
      };    
    };

    cout<<"Vamos en la linea i="<<i<<endl;
	
  };


  ZeroImag.close();
  testa.close();

  return 0;
  
}
