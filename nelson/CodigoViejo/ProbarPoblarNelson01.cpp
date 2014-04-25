//Veamos que podemos hacer con esto
//Tratemos de encontrar los Blind Spots Subplanckianos
//con esta chiva.


//No necesitas incluir nada que ya se haya incluido antes

//En esta version 04, la distribucion es radial. Asi sera facil comparar los
// ZEROS.


#include "RutinasNelson01.hpp"


using namespace std;

//Esto es pleonasmo, pero para que te acuerdes.

int main(){

  srand(23890);

  //Solo piensa en superposicion de ondas medio aleatoreas
  int muestreo=200;
  int nivel=300;
  mat DiracDeltas;

  DiracDeltas=PopulateNelson(nivel, muestreo);
  DiracDeltas.save("DiracNelson.dat", raw_ascii);
  


  return 0;
  

}
