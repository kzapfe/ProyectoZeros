//#include "../nelson/RutinasNelson02.hpp"


using namespace arma;


mat Pop(double energia, int muestreo){
  // Matrix de condiciones iniciales sobre el plano de simetria.

  double Kinetica, Potencial;
  mat result;
  int i=0;

    
  result=zeros<mat>(muestreo, 4);
  

   for(int i=0; i<muestreo;i++){
  
    Kinetica=energia*as_scalar(randu(1));
    Potencial=energia-Kinetica;
    rowvec auxiliar(4);
    double theta;
    //First Coin Toss: direccion del momento
    //Solo queremos x.p>0
    theta=-pi/2.0+pi*as_scalar(randu(1));
         
    //momentos 1 y 2
    auxiliar(1)=sqrt(2.0*Kinetica)*cos(theta);
    auxiliar(3)=sqrt(2.0*Kinetica)*sin(theta);

    auxiliar(0)=0.0000; //Plano de simetria
    
    auxiliar(2)=sqrt(Potencial);
    
    result.row(i)=auxiliar;

  };


  return result;

};

mat PopSmallRectangle(double energia, int muestreo){
  
  mat result;
  int i=0;

    
  result=zeros<mat>(muestreo, 4);


  
   for(int i=0; i<muestreo;i++){
     
  
     
     rowvec auxiliar(4);
     double cachito=0.001;
     //First Coin Toss: direccion del momento
     //Solo queremos x.p>0
     auxiliar(2)=cachito*as_scalar(randu(1));
      
     auxiliar(3)=cachito*as_scalar(randu(1));
       
     auxiliar(0)=0.0000; //Plano de simetria
     double energiarestante=energia-auxiliar(2)*auxiliar(2)
       -auxiliar(3)*auxiliar(3)/2.0;

     auxiliar(1)=sqrt(2.0*energiarestante);
     
     result.row(i)=auxiliar;
     
  };

  

  return result;

};



mat PopLineAtRandom(double energia, int muestreo){
  //Many mirrorsymetric orbits cross the x_q plane
  //With a perpendicular momentum
  mat result;
  int i=0;

    
  result=zeros<mat>(muestreo, 4);


  
   for(int i=0; i<muestreo;i++){
     
     rowvec auxiliar(4);
      double cachito=sqrt(energia);
      // double cachito=0.0200;
     
     
//First Coin Toss: direccion del momento
     //Solo queremos x.p>0
     //auxiliar(2)=-cachito+2.0*cachito*as_scalar(randu(1));[
      auxiliar(2)=0.13+cachito*as_scalar(randu(1));
      
     auxiliar(3)=0.00;
       
     auxiliar(0)=0.0000; //Plano de simetria
     double energiarestante=energia-auxiliar(2)*auxiliar(2)
       -auxiliar(3)*auxiliar(3)/2.0;

     auxiliar(1)=sqrt(2.0*energiarestante);
     
     result.row(i)=auxiliar;
     
  };

  

  return result;

};


mat PopLineInOrder(double energia, int muestreo){
  //Many mirrorsymetric orbits cross the x_q plane
  //With a perpendicular momentum
  mat result;
  int i=0;

    
  result=zeros<mat>(muestreo, 4);


  
   for(int i=0; i<muestreo;i++){
     
     rowvec auxiliar(4);
      double cachito=sqrt(energia);
      // double cachito=0.00000200;
     
     
//First Coin Toss: direccion del momento
     //Solo queremos x.p>0
     //auxiliar(2)=-cachito+2.0*cachito*as_scalar(randu(1));[
      auxiliar(2)=0.00+cachito*(double)i/(double)muestreo;
      
      auxiliar(3)=0.00;
      
      auxiliar(0)=0.0000; //Plano de simetria
      double energiarestante=energia-auxiliar(2)*auxiliar(2)
	-auxiliar(3)*auxiliar(3)/2.0;

      auxiliar(1)=sqrt(2.0*energiarestante);
      
      result.row(i)=auxiliar;
      
  };
   
  
  return result;

};
