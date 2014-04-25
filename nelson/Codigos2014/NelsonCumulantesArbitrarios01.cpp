/*
Expansion del operador de traslacion a la Fabricio.
en Cumulantes. 
 */

#include <armadillo>
#include <vector>
#include "simplectic01.hpp"
#include "ParametrosGlobales.hpp"
#include "RutinasNelson02.hpp"
#include "BinomialCoefficient01.hpp"
#include <gsl/gsl_sf_gamma.h>

using namespace std;

//Esto es pleonasmo, pero para que te acuerdes.

int main(){
  
  //inicializar el semillador
  srand(12389);

  const int talgrado=23;

  //Parametros Globales
 
  cout<<"la seccion tiene corte en mu (q ,p)= "<<endl;
  cout<<mu.q<<"\t"<<mu.p<<endl;


   
  std::ostringstream escupefuncion;
  std::string haz;
  escupefuncion<<"Nelson"<<"_"<<talgrado<<"grado"<<
    "_Weyl.dat"<<std::ends;
  haz=escupefuncion.str();
  const char *nombrefuncion=haz.c_str();
  ofstream Weyl;
  Weyl.open(nombrefuncion);

 
  mat DiracDeltas;
  DiracDeltas.load("CentrosWigner.dat");
  
    
  /*No pretendamos hacerlo tan general. Empecemos por el corte
  Que nos interesa: xhi tranformada de y (x2 en notacion Alfredo)
  mu transformada de x (x1 en Alf). mu=(0,0) -> 
  exp (-i/h (cuerdas wedge centros))=
  exp(-i/h xhi cuerda y) *1 
  Por ende solo nos interesan las dos ulimas columnas de DiracDeltas*/
  


  //para cada grado hacer los cumulantes necesarios
  //okey, te quieres ver muy pistola, not necesarry
  //no son tan grandes para nuestros propositos.
  
  //los grados de un polinomio empiecan en zero
  
  Mat<int> factor;
  mat cumulant;
  //no necesariamente tienen que corresponer a q y p
  vec auxq, auxp;
  
  auxq=DiracDeltas.col(2);
  auxp=DiracDeltas.col(3);
  

  factor.zeros(talgrado+1, talgrado+1);
  cumulant.zeros(talgrado+1, talgrado+1);

  for(int n=0; n<=talgrado; n++){
    //cout<<"vamos bien "<< j <<endl;
    for(int k=0; k<=n; k++){
      /*En este caso, por el corte escojido
	los elementos con potencias nones de auxp son negativos */
	
      factor(n,k)=ComputeBinomialCoefficient(n,k)*pow(-1,k);
	
      //cout<<factor(j,k)<<"\t";
      cumulant(n,k)=as_scalar(dot(pow(auxq, n-k),pow(auxp,k)))/(double)muestreo;
      cout<<cumulant(n,k)<<"\t";

    }
    cout<<endl;
  }

  cumulant.save("Cumulantes.dat", arma_ascii);
  factor.save("Binomiales.dat", arma_ascii);
 
  double weylreal=0.00;
  double weylimag=0.00;
  int enteroauxiliar;

  double xhi_q;
  double xhi_p;

  for(int i=-resolucion; i<resolucion; i++){
   
    for(int j=-resolucion; j<resolucion; j++){
      
       xhi_q=(double)i/(double)resolucion*extrem;
       xhi_p=(double)j/(double)resolucion*extrem;

       weylreal=0.00;
       weylimag=0.00;
    
  

       for(int n=0; n<=talgrado; n++){
	 for(int k=0; k<=n; k++){
	   enteroauxiliar=(n/2);
	   
	   if(n%2==0){
	     weylreal+=pow(xhi_p,n-k)*pow(xhi_q,k)*
	       pow(-1,enteroauxiliar)*factor(n,k)*
	       cumulant(n,k)/(gsl_sf_fact(n)*pow(hbarra,k));

	     
	   }else{
	     
	     weylimag+=pow(xhi_p,n-k)*pow(xhi_q,k)*
	       pow(-1,enteroauxiliar+1)*factor(n,k)*
	       cumulant(n,k)/(gsl_sf_fact(k)*pow(hbarra,k));	      
	  
	   }//par o non
	  

	 } //sobre los elementos de la suma
	 
       } //sobre el grado de la aproximacion
       
  
       
       Weyl<<xhi_q<<"\t"<<xhi_p<<"\t"<<
	 weylreal<<"\t"<<weylimag<<endl;
      

    };

    
    
    Weyl<<endl;

    // cout<<"Vamos en el i="<<i<<endl;
  };


  Weyl.close();
  return 0;

}
