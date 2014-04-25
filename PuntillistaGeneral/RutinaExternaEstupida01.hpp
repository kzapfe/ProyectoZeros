/*Las secciones en centros y cuerdas */

using namespace std;
using namespace arma;

void TestRecibeVector(vector<simplectic>& Xint, vector<simplectic>& Yint, 
		      string nombrearchivo){
  

  
  ofstream Archivo;
  //Funcion chingona que convierte un string a c_str que funciona para open.
  Archivo.open(nombrearchivo.c_str()); 
  Archivo<<"#Test" <<endl;
  
  
  int maximumpoints;
  maximumpoints=Xint.size();
  
  
    
  for(int i=0; i<maximumpoints; i ++){
    
    Archivo<<Xint[i].q<<"\t"<<Yint[i].p<<endl;
      
      
      }
	
    
  
  Archivo.close();
  
  //Aqui termina su trabajo el hilo 0
  
  
  return;

};





