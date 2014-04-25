//Veamos que podemos hacer con esto
//Tratemos de encontrar los Blind Spots Subplanckianos
//con esta chiva.


//No necesitas incluir nada que ya se haya incluido antes

//En esta version 04, la distribucion es radial. Asi sera facil comparar los
// ZEROS.


#include "RutinasNelson01.hpp"
#include "intervalito02.hpp"

using namespace std;

//Esto es pleonasmo, pero para que te acuerdes.

int main(){
  
  //inicializar el semillador
  srand(23890);

  //Parametros Globales
  int muestreo=200;
  int nivel=300;

  const int resolucion1=3; //La resolucion es una propiedad de la rutina.
  const int resolucion2=5; //La resolucion es una propiedad de la rutina.
  const int resolucion3=5; //La resolucion es una propiedad de la rutina.
  

  //Abrimos los archivos necesarios
  std::ostringstream escupezeros;
  escupezeros<<muestreo<<"_NelsonZerosReales.dat"<<std::ends;
  std::string stringzeros;
  stringzeros=escupezeros.str();
  const char *nombrezerosreal=stringzeros.c_str();
  ofstream zerolinereal;
  zerolinereal.open(nombrezerosreal);

 
  escupezeros<<muestreo<<"_NelsonZerosImag.dat"<<std::ends; 
  stringzeros=escupezeros.str();
  const char *nombrezerosimag=stringzeros.c_str();
  
  ofstream zerolineimag;
  zerolineimag.open(nombrezerosimag);

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
  
  //----------------------------------------------------
  
  //********HERE BEGINS THE PARTY !!!******************
  //----------------------------------------------------




  for(int m=0; m<3;m++){ 
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
	      int bigrhorealsign=0;
	      int smallrhorealsign=0;
	      int bigrhoimagsign=0;
	      int smallrhoimagsign=0;
	      
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
	     
	     


	     smallrhorealsign=signum(weylreal);     
	     smallrhoimagsign=signum(weylimag);     

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
	  


	     bigrhorealsign=signum(weylreal);
	     bigrhoimagsign=signum(weylimag);
  
	     //**************************************************	  
	     
	     result=(bigrhorealsign!=smallrhorealsign);
	    
	     subintervalos[pieza].cambiareal=result;

	     result=(bigrhorealsign!=smallrhoimagsign);

	     subintervalos[pieza].cambiaimag=result;

	     //FIN DA MERDE



	    
	      if(subintervalos[pieza].cambia){
		auxinterval=subintervalos[pieza];
		auxinterval.cambia=true;
		

		pieza=divisiones+1;
		iteracion++;
		
	      }else{
		pieza++;
		
		if(pieza==divisiones){		 
		  iteracion++;
		};
	      }; 
	      
	  };
	    
	};
	
	zerolinereal<<auxinterval.menor<<"\t"<<theta1<<"\t"
		<<theta2<<"\t"
		<<theta3<<"\t"
		<<iteracion<<"\t"
		<<auxinterval.cambia
		<<endl;
	    	  
      };

      };
      zerolinereal<<endl;
      zerolinereal<<endl;
  };
  
  };

  
  zerolinereal.close();
 

  return 0;
  

}
