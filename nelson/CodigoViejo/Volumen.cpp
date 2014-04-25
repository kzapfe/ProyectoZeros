//Obtener el volumen de la camada de energia.

#include "RutinasNelson02.hpp"
#include "ParametrosGlobales.hpp"

using namespace std;



int main(){  
  //inicializar el semillador
  srand(12389);  

  
  
  ofstream testa;
  testa.open("Datosparachecar.dat");
  
  
  //----------------------------------------------------
  
  //Poblar la camada de energia
  mat DiracDeltas;
  DiracDeltas=PopulateNelson(nivel, muestreo);
  
  double distprom=0.000;
  double distanciaacumulada=0.00;
  double estimatedvol=0.000;

  for(int i=0; i<muestreo; i++){
    for(int j=0; j<muestreo; j++){
      //La idea trivial NO funciona
      //Si el conjunto no es convexo.
      //Piensa dos pelotas separadas por una distancia
      //gigante, y chan chan.
      distanciaacumulada+=
	distancia(DiracDeltas.row(i),DiracDeltas.row(j));
    };
    cout<<i<<endl;
  };
  

  distprom=distanciaacumulada/(2.0*muestreo);
  estimatedvol=pow(distprom, 4.0);

  cout<<"La distancia promedio para cada punto es "
      <<distprom<<endl;
  cout<<"El Volumen asi estimado es "<<estimatedvol<<endl;

  return 0;
  

}

      
    
