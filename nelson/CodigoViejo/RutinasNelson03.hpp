//esta version 03 utiliza Una distribucion que cae como una maxwelliana
// fuera de la camada de energia. 





using namespace arma;
using namespace std;

int signum(double x){
  int sign;
  if(x>0){
    sign=1;
  }else if(x<0){
    sign=-1;
  }else{
    sign=0;
      };
  
  return sign;  
};


double NelsonEnergy(rowvec x){

  double result;
  double q1,q2,p1,p2;
  
  q1=x(0);
  p1=x(1);
  q2=x(2);
  p2=x(3);
  
   result=(p1*p1+p2*p2)/2.0+
    q1*q1/20.0
    +(q2-q1*q1/2)*(q2-q1*q1/2);
   
  return result;

};


double ValorExactoq2(double q1, double  Potencial){
  
  //Aqui tambien tiramos un dado
  double result;
  double aux=0.000;
  aux=as_scalar(randn(1));
  if(aux>0.16700){
    result=(q1*q1+sqrt(4.0*Potencial-q1*q1/5.0))/2.0;
  }else{
    result=(q1*q1-sqrt(4.0*Potencial-q1*q1/5.0))/2.0;
  };
  
  return result;

};

double DistLineal(){

  //Producimos un x aleatoreo
  //con probabilidad diferencial 1-x
  //a partir de un y homogeneamente generado dp/dy=1
  double y,x;
  y=as_scalar(randu(1));
  x=1-sqrt(1-y);
  return x;

};

//Use your head Karel
//To populate a BALL you dont care if it is symplectic!
//Actually, in your approach,you SCREW simplecticity!!
rowvec CreateRandomPointinBall(double nivel=0, int dimension=2){
  /*Montecarleando la capa de energia numero nivel
  //Mete uno con distribucion Gaussiana en la energia
  Y distribucion Uniforme en el angulo */

  //Marsaglia algoritm
  rowvec x(dimension);
  double radio, angulo;
  double centro; 
  double delta;
  centro=sqrt(hbarra * nivel);
  delta=sqrt(hbarra * (nivel+1))-sqrt(hbarra * nivel);

  x.randn();
  x=x/norm(x,2);
  // aqui ya tenemos un punto en la esfera de radio 1

  radio=as_scalar(randu(1))*delta;
  radio=pow(radio+centro, 1.0/(double)dimension);
  radio=radio*pow(delta+centro, (double)(dimension-1)/(double)dimension);
  
  x=x*radio;
  
  return x;
  
};


mat Populate4BallShell(double nivel=0, int muestreo=1){
  //Crea un array del tamano de muestreo
  // poblando gausianamente cerca de la capa de energia que
  // nos interesa, dada por nivel.
  
  mat points(muestreo, 4);
  

  for(int i=0; i<muestreo; i++) {
    points.row(i)= CreateRandomPointinBall(nivel, 4);
  };

  return points;
    
};


mat PopulateNelson(double nivel=0, int muestreo=1){
  //Usa la cabezota, guey
  //para cada valor de E, q2 y |p| pueden ser exactamente determinados.
  //entonces realmente lo unico que necesitamos
  //es escoger aleatoreamente TRES de las cuatro, y variarle una
  //unidad a E.

  double particion;
  double Energia, Potencial, Kinetica;
  mat result;
  double Ham;
  int i=0;
  double theta;
  double xaux,yaux;

  //Tu unidad de energia mas pequenha es 0.05*hbarra (freq. menor)
    
  double q1maxvalue;
  double q2max,q2min;

  double sigma;
  double centermu;
  centermu=nivel;
  sigma=hbarra*0.1/4.0;
  
  result=zeros<mat>(muestreo, 4);

  rowvec auxiliar=zeros<rowvec>(4);

  for(int i=0; i<muestreo;i++){
    

    //First Coin Toss: lugar en el gordito entre nivel y  nivel+1 
   
    Energia=as_scalar(randn(1))*sigma+centermu;
    if(Energia<0.00)Energia=-Energia;    


    //la distribucion estaba fea
    //probemos lo siguiente: q_1 va a estar aleatoreamente
    //entre las posiciones maximas permitidas.
    q1maxvalue=sqrt(20.0*Energia); 
    //second coin toss
    //THIS MAKES IT TO ACUMULATE SLIGHTLY ON THE HORNS
    //auxiliar(0)=as_scalar(randu(1))*2.0*q1maxvalue-q1maxvalue;
    //    cout<< "q1 es= "<<auxiliar(0)<<endl;

    //Vamos a tratar de darle un peso lineal

    xaux=DistLineal();
    yaux=as_scalar(randu(1))-0.5;
    xaux=xaux*signum(yaux);

    auxiliar(0)=xaux*q1maxvalue;

    //q2 va a estar entre su valor maximo y minimo permitido

    q2min=(auxiliar(0)*auxiliar(0)
	   -sqrt(4.0*Energia-auxiliar(0)*auxiliar(0)/5.0))/2.0;
    q2max=(auxiliar(0)*auxiliar(0)
	   +sqrt(4.0*Energia-auxiliar(0)*auxiliar(0)/5.0))/2.0;

    //third coin toss
    auxiliar(2)=as_scalar(randu(1))*(q2max-q2min)+q2min;
    //cout<< "q2 es= "<<auxiliar(2)<<endl;


    //Ver cuanta energia nos queda libre
    Potencial=auxiliar(0)*auxiliar(0)/20.0+
      (auxiliar(2)-auxiliar(0)*auxiliar(0)/2.0)*
      (auxiliar(2)-auxiliar(0)*auxiliar(0)/2.0);
    
    //cout<< "Potencial es= "<<Potencial<<endl;

    Kinetica=Energia-Potencial;


    //cuarto Coin Toss: direccion del momento
    theta=2.00*pi*as_scalar(randu(1));
        
        //momentos 1 y 2
    auxiliar(1)=sqrt(2.0*Kinetica)*cos(theta);
    
    auxiliar(3)=sqrt(2.0*Kinetica)*sin(theta);

    result.row(i)=auxiliar;
   
  };
      
  //  cout<<result<<endl;

  return result;

};



mat PopulateNelsonwithGauss(double nivel=0, int muestreo=1){
  //Usa la cabezota, guey
  //para cada valor de E, q2 y |p| pueden ser exactamente determinados.
  //entonces realmente lo unico que necesitamos
  //es escoger aleatoreamente TRES de las cuatro
  //vamos a tomar en cuenta la simetria especular.
  
  //simetrizamos el numero de entradas.
  if(muestreo%4 !=0) muestreo=muestreo-(muestreo%4);

  cout<<"malico"<<endl;
  cout<<muestreo<<endl;

  double particion;
  double Energia, Potencial, Kinetica;
  mat result;
  double Ham;
  int i=0;
  double theta;
  double xaux,yaux;

  
  
  //parametros de la distribucion
  double sigma;
  double centermu;

  //Tu unidad de energia mas pequenha es 0.1*hbarra (freq. menor)
  
  centermu=nivel;
  sigma=hbarra*0.1;
  
  double q1maxvalue;
  double q2max,q2min;

  
  result=zeros<mat>(muestreo, 4);

  rowvec auxiliar=zeros<rowvec>(4);

  for(int i=0; i<muestreo/4;i++){
  

    //First Coin Toss: lugar en el gordito entre nivel+1 y  nivel-1 
    Energia=as_scalar(randn(1))*sigma+centermu;
    if(Energia<0.00)Energia=-Energia;    

    //la distribucion estaba fea
    //probemos lo siguiente: q_1 va a estar aleatoreamente
    //entre las posiciones maximas permitidas.
    q1maxvalue=sqrt(20.0*Energia); 
    //second coin toss
    //THIS MAKES IT TO ACUMULATE SLIGHTLY ON THE HORNS
    //auxiliar(0)=as_scalar(randu(1))*2.0*q1maxvalue-q1maxvalue;
    //    cout<< "q1 es= "<<auxiliar(0)<<endl;

    //Vamos a tratar de darle un peso lineal

    xaux=DistLineal();
    yaux=as_scalar(randu(1))-0.5;
    xaux=xaux*signum(yaux);

    auxiliar(0)=xaux*q1maxvalue;

    //q2 va a estar entre su valor maximo y minimo permitido

    q2min=(auxiliar(0)*auxiliar(0)
	   -sqrt(4.0*Energia-auxiliar(0)*auxiliar(0)/5.0))/2.0;
    q2max=(auxiliar(0)*auxiliar(0)
	   +sqrt(4.0*Energia-auxiliar(0)*auxiliar(0)/5.0))/2.0;

    //third coin toss
    auxiliar(2)=as_scalar(randu(1))*(q2max-q2min)+q2min;
    //cout<< "q2 es= "<<auxiliar(2)<<endl;


    //Ver cuanta energia nos queda libre
    Potencial=auxiliar(0)*auxiliar(0)/20.0+
      (auxiliar(2)-auxiliar(0)*auxiliar(0)/2.0)*
      (auxiliar(2)-auxiliar(0)*auxiliar(0)/2.0);
    
    //cout<< "Potencial es= "<<Potencial<<endl;

    Kinetica=Energia-Potencial;


    //cuarto Coin Toss: direccion del momento
    theta=2.00*pi*as_scalar(randu(1));
        
        //momentos 1 y 2
    auxiliar(1)=sqrt(2.0*Kinetica)*cos(theta);
    
    auxiliar(3)=sqrt(2.0*Kinetica)*sin(theta);

    result.row(i)=auxiliar;
    //reflejo q_x--> -q_x
    auxiliar(0)=-auxiliar(0);
    result.row(i+muestreo/4)=auxiliar;
    //reflexo p_y --> -p_y y la anterior
    auxiliar(3)=-auxiliar(3);
    result.row(i+muestreo/2)=auxiliar;
    //reflexo p_y --> -p_y nadamas
    auxiliar(0)=-auxiliar(0);
    result.row(i+3*muestreo/4)=auxiliar;
    
    
 
  };
    

  //  cout<<result<<endl;

  return result;

};



double WeylReal(int muestreo, simplectic &cuerda1, simplectic &cuerda2,
		simplectic * centro1, simplectic * centro2){

  double result=0.00;  
  for(int k=0; k<muestreo; k++){
    //Evaluando extremo inferior del intervalo
    
    result+=cos(-(centro1[k].simplecticproduct(cuerda1)+
			  centro2[k].simplecticproduct(cuerda2))/hbarra);
    
  };

  result=result/(hbarra*pi*muestreo);  
  return result;

}; 


double WeylImag(int muestreo, simplectic &cuerda1, simplectic &cuerda2,
		simplectic * centro1, simplectic * centro2){

  double result=0.00;  
  for(int k=0; k<muestreo; k++){
    //Evaluando extremo inferior del intervalo
    
    result+=sin(-(centro1[k].simplecticproduct(cuerda1)+
			  centro2[k].simplecticproduct(cuerda2))/hbarra);
    
  };

  result=result/(hbarra*pi*muestreo);  
  return result;

}; 



double distancia(rowvec x, rowvec y){
  //distancia euclides entre dos rowvec
  double result=0;

  result=norm((x-y),2);
  
  return result;

};


mat PopulateCircle(double nivel=0, int muestreo=1){
  //Crea un array del tamano de muestreo
  // poblando la cicatriz correspondiente a la primera 
  // Orbita periodica del potencial de Nelson.
  
  mat points(muestreo, 4);
  double FabrizioEpsilon=0.0025;
  double energiaconruido=nivel
    +FabrizioEpsilon*as_scalar(randu(1));


  for(int i=0; i<muestreo; i++) {

    rowvec auxiliar(4);
    double theta;
  
    theta=2.00*pi*as_scalar(randu(1));
            
    //Componentes y
    auxiliar(2)=sqrt(energiaconruido)*cos(theta);
    
    auxiliar(3)=sqrt(energiaconruido*2.00)*sin(theta);
  

    //componentes x
    auxiliar(0)=0.00;

    auxiliar(1)=0.000;
    
    //    double potaux=auxiliar(0)*auxiliar(0)/20.0+
    // (auxiliar(2)-auxiliar(0)*auxiliar(0)/2.0)*
    // (auxiliar(2)-auxiliar(0)*auxiliar(0)/2.0);
    

    points.row(i)=auxiliar;
    
  };

    
 

  
  return points;
    
};

mat PopulatefromDrive(int muestreo, string archivo){
  //Tratemos de escoger un buen muestro a
  //partir de datos en el disco duro
  mat data;
  data.load(archivo, auto_detect);
  mat temp=shuffle(data);
  mat result;

  cout<<"pito"<<endl;
  result=temp.submat(0,0,muestreo-1,3);
  cout<<"pito2"<<endl;

  return result;

};


