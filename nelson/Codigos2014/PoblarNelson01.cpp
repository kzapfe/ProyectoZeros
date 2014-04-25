/*
Solamente poblar la camada de energia y guardar los Centros Wigner
This will be used again for the Cat States 
*/

#include <armadillo>
#include "simplectic01.hpp"
#include "ParametrosGlobales.hpp"
#include "RutinasNelson02.hpp"


using namespace std;

//Esto es pleonasmo, pero para que te acuerdes.

int main(){
  
  //inicializar el semillador
  srand(12389);
  int muestraparacomparar=1400;
 
  mat DiracDeltas;
  DiracDeltas=PopulateNelson(nivel, muestraparacomparar);
  DiracDeltas.save("CentrosNelson.dat", arma_ascii);
    

  return 0;

}
