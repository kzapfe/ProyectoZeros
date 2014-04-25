//Veamos que podemos hacer con esto
//Tratemos de encontrar los Blind Spots Subplanckianos
//con esta chiva.


//No necesitas incluir nada que ya se haya incluido antes

//En esta version 04, la distribucion es radial. Asi sera facil comparar los
// ZEROS.


#include "PuntosCuerdasEsfera01.hpp"


using namespace std;

//Esto es pleonasmo, pero para que te acuerdes.

int main(){

  srand(23890);

  //Solo piensa en superposicion de ondas medio aleatoreas
  int muestreo;
  int nivel=300;

  /*
  muestreo=500;
  CalculaWeyl(muestreo, nivel);
  */

    
  muestreo=20;
  CalculaWeyl(muestreo, nivel);
  cout<<"Sale baja Res"<<endl;
  muestreo=200;
  CalculaWeyl(muestreo, nivel);
  cout<<"Sale media Res"<<endl;
  /*  muestreo=2000;
  CalculaWeyl(muestreo, nivel);
  cout<<"Sale alta Res"<<endl;
  muestreo=20000;
  CalculaWeyl(muestreo, nivel);
  cout<<"Sale altisima Res"<<endl;

  */

  return 0;



}
