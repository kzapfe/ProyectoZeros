/*Vamos a buscar las orbitas de mas bajo
  periodo usando dinamica clasica */

#include <armadillo>
#include "../nelson/simplectic01.hpp"
#include "../nelson/ParametrosGlobales.hpp"
#include "RutinasAuxiliares01.hpp" 
#include "DerivadasHamilton01.hpp" 


using namespace std;

//using std::setprecision;



int main(){
  
  simplectic x,y;
  double energia=0.821;
  double epsilon=0.0001;
  int Time=5000000;
  int kmax=2;

  ofstream Poincare;
  ofstream Orbits;
  int numerodeorbitas=40000;
  mat condini=zeros<mat>(numerodeorbitas,4);

  streamsize prec=Poincare.precision(12);


  Poincare.open("PoincareSample.dat");
  Orbits.open("Orbitas02.dat");

  condini=PopLineInOrder(energia, numerodeorbitas);
  
  condini.save("Condicionesiniciais.dat", raw_ascii);
  
  for(int i=0; i<numerodeorbitas; i++){
    x.q=condini(i,0);
    x.p=condini(i,1);
    y.q=condini(i,2);
    y.p=condini(i,3);

    int kontamaps=0;
    
    Orbits<<x.q<<"\t"<<x.p<<"\t"
	  <<y.q<<"\t"<<y.p<<"\t"
	  <<0<<"\t"<<i<<endl;
    

    Poincare<<x.q<<"\t"<<x.p<<"\t"
	    <<y.q<<"\t"<<y.p<<"\t"
	    <<0<<"\t"<<kontamaps
	    <<"\t"<<i<<endl;	  
    
    for(int t=1; t<Time; t++){
      double xqaux;
      double yqaux;

      xqaux=x.q;
      yqaux=y.q;
      //Vuelo libre infinitesimal
      x.q+=epsilon*x.p;
      y.q+=epsilon*y.p;
	
      //PatadaInfinitesimal
      x.p-=epsilon*dH_dxq(x,y);
      y.p-=epsilon*dH_dyq(x,y);
      
      if((t%1000==0)&&(i%40==0)){
	Orbits<<x.q<<"\t"<<x.p<<"\t"
	      <<y.q<<"\t"<<y.p<<"\t"
	      <<t*epsilon<<"\t"<<i<<endl;
      };
      
      if((xqaux<0.00000)&&(x.q>0.00000)&&(t>5)){
	kontamaps++;
	double xqinterpol, yqinterpol;
	xqinterpol=0.0000;
	yqinterpol=yqaux+y.p*(-xqaux/x.p);
      
	Poincare<<xqinterpol<<"\t"<<x.p<<"\t"
		<<yqinterpol<<"\t"<<y.p<<"\t"
		<<t*epsilon<<"\t"<<kontamaps
		<<"\t"<<i<<endl;
      };
      
      if(kontamaps==kmax)t=Time;
      

    };
	
    Orbits<<endl;
    Poincare<<endl;

  };


  return 0;


}
