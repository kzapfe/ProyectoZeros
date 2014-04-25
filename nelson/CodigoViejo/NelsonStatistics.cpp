/*
Es te programa dibuja la funcion de Cuerdas para 
planos mu=cte, en subdominio xhi. 

Version 02: Vamos a poner este desmadre en orden


 */
#include <armadillo>
#include "simplectic01.hpp"
#include "ParametrosGlobales.hpp"
#include "RutinasNelson02.hpp"
#include "ChordAprox01.hpp"


using namespace std;

//Esto es pleonasmo, pero para que te acuerdes.

int main(){
  
  //inicializar el semillador
  srand(12389);

  //Parametros Globales
 
  cout<<"la seccion tiene q y p: "<<endl;
  cout<<mu.q<<"\t"<<mu.p<<endl;
  //Abrimos los archivos necesarios

 

  std::ostringstream escupemedias;
  std::string hazstring;
  escupemedias<<"ExactoN"<<"_"<<"0821"<<
    "_CumulantesWigner0-0.dat"<<std::ends;
  hazstring=escupemedias.str();
  const char *nombremedwig=hazstring.c_str();
  ofstream MediasWigner;
  MediasWigner.open(nombremedwig);



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

  vec nicemean;
  nicemean=conv_to<colvec>::from(meanvalues);
  
  MediasWigner<<"Primer momento x_q: "<< nicemean(0)<<endl;
  MediasWigner<<"Primer momento x_p: "<< nicemean(1)<<endl;
  MediasWigner<<"Primer momento y_q: "<< nicemean(2)<<endl;
  MediasWigner<<"Primer momento y_p: "<< nicemean(3)<<endl;

  MediasWigner<<endl;
  

  meanvalues=var(DiracDeltas, 1);

   
  MediasWigner<<"Segundo momento x_q: "<< meanvalues(0)<<endl;
  MediasWigner<<"Segundo momento x_p: "<< meanvalues(1)<<endl;
  MediasWigner<<"Segundo momento y_q: "<< meanvalues(2)<<endl;
  MediasWigner<<"Segundo momento y_p: "<< meanvalues(3)<<endl;
  
  MediasWigner<<endl;
  


  mat MatrixKovariancia=zeros<mat>(4,4);
  
  //This is easy
  for(int i=0; i<4;i++){
    for(int j=0; j<4;j++){
      MatrixKovariancia(i,j)=
	dot(DiracDeltas.col(i),DiracDeltas.col(j))/(double)muestreo;
    };
  };
  
  MediasWigner<<" Esta el la matrix de Kovariancia K: "<<endl;
  MediasWigner<<MatrixKovariancia<<endl;

  MediasWigner<<endl;

  //Up to the second non central cumulant.

    
  cube KurtosisCube=zeros<cube>(4,4,4);
  //Imagina el vector (xhi_q,xhi_p, mu_q,mu_p)
  //Usa ese orden par numerar los coeficientes
  //Exampligratia: el coeficiente de xhi_q *mu_q *mu_p es
  // KurtosisCube(0,2,3) Armadillo cuenta desde cero, recuerda.
  //Dadas simetrias, esto esta de hueva pero sale rapido
  
  //Primero, los que multiplican a las potencias puras
  //aguado con los signos
  KurtosisCube(0,0,0)=
    -dot(DiracDeltas.col(1)%DiracDeltas.col(1),DiracDeltas.col(1))
    /(double)muestreo;
  KurtosisCube(1,1,1)=
    dot(DiracDeltas.col(0)%DiracDeltas.col(0),DiracDeltas.col(0))
    /(double)muestreo;
  KurtosisCube(2,2,2)=
    -dot(DiracDeltas.col(3)%DiracDeltas.col(3),DiracDeltas.col(3))
    /(double)muestreo;
  KurtosisCube(3,3,3)=
    dot(DiracDeltas.col(2)%DiracDeltas.col(2),DiracDeltas.col(2))
    /(double)muestreo;
  
  //Ahora los que tienen uno de dos y uno de uno
  //Del mismo par simplectico
  KurtosisCube(0,0,1)=
    3.0*dot(DiracDeltas.col(1)%DiracDeltas.col(1),DiracDeltas.col(0))
    /(double)muestreo;
  KurtosisCube(0,1,1)=
    -3.0*dot(DiracDeltas.col(0)%DiracDeltas.col(0),DiracDeltas.col(1))
    /(double)muestreo;
  KurtosisCube(2,2,3)=
    3.0*dot(DiracDeltas.col(3)%DiracDeltas.col(3),DiracDeltas.col(2))
    /(double)muestreo;
  KurtosisCube(3,3,2)=
    -3.0*dot(DiracDeltas.col(2)%DiracDeltas.col(2),DiracDeltas.col(3))
    /(double)muestreo;
  
  //Ahora los de dos de uno y uno de otro de distintos pares simplecticos

  // qi qi qj
  KurtosisCube(0,0,2)=
    -3.0*dot(DiracDeltas.col(1)%DiracDeltas.col(1),DiracDeltas.col(3))
    /(double)muestreo; 
  KurtosisCube(2,2,0)=
    -3.0*dot(DiracDeltas.col(3)%DiracDeltas.col(3),DiracDeltas.col(1))
    /(double)muestreo;
 
  // q q p
  KurtosisCube(0,0,3)=
    3.0*dot(DiracDeltas.col(1)%DiracDeltas.col(1),DiracDeltas.col(2))
    /(double)muestreo;
  KurtosisCube(2,2,1)=
    3.0*dot(DiracDeltas.col(3)%DiracDeltas.col(3),DiracDeltas.col(0))
    /(double)muestreo;
  
  // p p p
   KurtosisCube(1,1,3)=
    3.0*dot(DiracDeltas.col(0)%DiracDeltas.col(0),DiracDeltas.col(2))
    /(double)muestreo;
   KurtosisCube(3,3,1)=
     3.0*dot(DiracDeltas.col(2)%DiracDeltas.col(2),DiracDeltas.col(0))
     /(double)muestreo;
   
   // p p q
   KurtosisCube(1,1,2)=
    -3.0*dot(DiracDeltas.col(0)%DiracDeltas.col(0),DiracDeltas.col(3))
    /(double)muestreo;
   KurtosisCube(3,3,0)=
     -3.0*dot(DiracDeltas.col(2)%DiracDeltas.col(2),DiracDeltas.col(1))
     /(double)muestreo;
  
   
   //Los que tienen una de cada una
   
   // q1 q2 p1
   KurtosisCube(1,1,3)=
     -6.0*dot(DiracDeltas.col(0)%DiracDeltas.col(0),DiracDeltas.col(2))
     /(double)muestreo;
   KurtosisCube(3,3,1)=
     -6.0*dot(DiracDeltas.col(2)%DiracDeltas.col(2),DiracDeltas.col(0))
     /(double)muestreo;
  
   // q1 p1 q2
   KurtosisCube(0,1,2)=
     6.0*dot(DiracDeltas.col(1)%DiracDeltas.col(0),DiracDeltas.col(3))
     /(double)muestreo;
   KurtosisCube(2,3,0)=
     6.0*dot(DiracDeltas.col(3)%DiracDeltas.col(2),DiracDeltas.col(1))
     /(double)muestreo;
  

   MediasWigner<< "El tensor de chuequez "<<endl;
   MediasWigner<<KurtosisCube<<endl;

   // cout<<KurtosisCube<<endl;

    
   vec chordvector;
   chordvector=zeros<vec>(4);
   chordvector(1)=0.5;

   double test;
   test=XiUno(nicemean, chordvector);
   cout<< "eyes punto median = "<<test<<endl;
   double test2=XiDos(MatrixKovariancia, chordvector);
   cout<< "eyes K eyes = "<<test2<<endl;
   double test3=XiTres(KurtosisCube, chordvector);
   cout<< "esyes Kurt eyes eyes = "<<test3<<endl;
  
   cout<<" Ahora hagamos la funcion con estos aproximantes!"<<endl;


   

  std::ostringstream escupefuncion;
  std::string haz;
  escupefuncion<<"ExactoN"<<"_"<<"0821"<<
    "_WeylAprox3grado-0-0.dat"<<std::ends;
  haz=escupefuncion.str();
  const char *nombrefuncion=haz.c_str();
  ofstream Weyl;
  Weyl.open(nombrefuncion);

  double weylreal=0.00;
  double weylimag=0.00;
  

  for(int i=-resolucion; i<resolucion; i++){
    for(int j=-resolucion; j<resolucion; j++){
      double xhi_q=(double)i/(double)resolucion*extrem;
      double xhi_p=(double)j/(double)resolucion*extrem;
      chordvector(0)=mu.q;
      chordvector(1)=mu.p;
      chordvector(2)=xhi_q;
      chordvector(3)=xhi_p;

      weylreal=1.00-XiDos(MatrixKovariancia, chordvector)/(2.0*hbarra*hbarra);

      weylimag=-XiUno(nicemean, chordvector)/hbarra+
	XiTres(KurtosisCube, chordvector)/(6.0*hbarra*hbarra*hbarra);
      
      weylreal=weylreal;
      weylimag=weylimag;

      Weyl<<xhi_q<<"\t"<<xhi_p<<
	"\t"<<weylreal<<"\t"<<weylimag<<endl;

    };
    
    Weyl<<endl;

    // cout<<"Vamos en el i="<<i<<endl;
  };

  MediasWigner.close();  
  Weyl.close();
  return 0;

}
