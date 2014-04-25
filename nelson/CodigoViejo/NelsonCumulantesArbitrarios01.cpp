/*
Expansion del operador de traslacion a la Fabricio.
en Cumulantes. 
 */

#include <armadillo>
#include <vector>
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
 
  cout<<"la seccion tiene corte en mu (q ,p)= "<<endl;
  cout<<mu.q<<"\t"<<mu.p<<endl;
  //Abrimos los archivos necesarios

  std::ostringstream escupemedias;
  std::string hazstring;
  escupemedias<<"Nelson"<<"_"<<"0821"<<
    "_CumulantesWigner0-0.dat"<<std::ends;
  hazstring=escupemedias.str();
  const char *nombremedwig=hazstring.c_str();
  ofstream MediasWigner;
  MediasWigner.open(nombremedwig);

   
  std::ostringstream escupefuncion;
  std::string haz;
  escupefuncion<<"Nelson"<<"_"<<"0821"<<
    "_WeylAprox3grado-0-0-newmetod.dat"<<std::ends;
  haz=escupefuncion.str();
  const char *nombrefuncion=haz.c_str();
  ofstream Weyl;
  Weyl.open(nombrefuncion);

 
  //----------------------------------------------------
  //Poblar la camada de energia
  mat DiracDeltas;
  DiracDeltas=PopulateNelson(nivel, muestreo);
  
  
  
  /*No pretendamos hacerlo tan general. Empecemos por el corte
  Que nos interesa: xhi tranformada de y (x2 en notacion Alfredo)
  mu transformada de x (x1 en Alf). mu=(0,0) -> 
  exp (-i/h (cuerdas wedge centros))=
  exp(-i/h xhi cuerda y) *1 */
  
  const int talgrado=3;

  double weylreal=0.00;
  double weylimag=0.00;
  

  for(int i=-resolucion; i<resolucion; i++){
    for(int j=-resolucion; j<resolucion; j++){
      
      double xhi_q=(double)i/(double)resolucion*extrem;
      double xhi_p=(double)j/(double)resolucion*extrem;
      
      xhi.simplecticproduct(y[k])

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
