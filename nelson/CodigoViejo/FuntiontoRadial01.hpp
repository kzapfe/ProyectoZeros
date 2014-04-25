/* Funcion pecueca que toma una matriz de datos
   representando una funcion compleja
   toma las dos primeras entras como la chingadera compleja
   y las otras dos como el valor
   y transforma el numerito a polares
   y las ordena */

#include <armadillo>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

using namespace std;
using namespace arma;


mat Radializador(mat FuncionComplex){
  mat matrizin;
  mat matrizout;
  int renglones=matrizin.n_rows;

  matrizout=matrizin;

  for(int i=0; i<renglones;i++){
    gsl_complex zhi;
    zhi=gsl_complex_rect(matrizin(i,0),matrizin(i,1));
    matrizout(i,0)=gsl_complex_abs(zhi);
    matrizout(i,1)=gsl_complex_arg(zhi);
  };
  
  return matrizout;


};
