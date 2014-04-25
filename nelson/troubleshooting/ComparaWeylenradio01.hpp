/*Esta rutina Solamente hace lo siguente:
  Fijas las variables angulares
  compara el valor de la funcion de Weyl
  en los dos extremos de un intervalo
  RADIAL.
 */



bool ComparaWeylRadial2D
(intervalito rad, double theta1, double theta2, double theta3, 
 simplectic *x, simplectic *y, int muestreo){
  

  bool result;
  
  double rho;
  simplectic xhi, mu;
  int bigrhoWeylSign=0;
  int smallrhoWeylSign=0;
  
  double weylreal=0.0;
  double weylimag=0.0;
   

  //**************************************************
  //Get the diabolical Weyl Function
  //for this mu and xhi
  rho=rad.menor;
  
  mu.q= cos(theta1)*rho;
  mu.p= sin(theta1)*cos(theta2)*rho;
  xhi.q=sin(theta1)*sin(theta2)*cos(theta3)*rho;
  xhi.p=sin(theta1)*sin(theta2)*sin(theta3)*rho;
  
  //Temporarly this merde here, cause if not, doesnt work.

  for(int k=0; k<muestreo; k++){
    //Aqui hacemos la transformada de Fourier
	    // Que nos da la Funcion de Weyl
    
    weylreal+=cos((x[k].simplecticproduct(mu)+
		   y[k].simplecticproduct(xhi))/hbarra);
    
    weylimag+=sin((x[k].simplecticproduct(mu)+
		   y[k].simplecticproduct(xhi))/hbarra);
    
  };
  
  
  
  weylreal=weylreal/(hbarra*pi*muestreo);
  weylimag=weylimag/(hbarra*pi*muestreo);
  
  //Weyl=WeylFunction2D(muestreo, x, y, mu, xhi);
  
  //weylreal=GSL_REAL(Weyl);
  //weylimag=GSL_IMAG(Weyl);



  smallrhoWeylSign=signum(weylreal);
  

  //**************************************************
  
  //**************************************************
  //Get the diabolical Weyl Function
  //for this mu and xhi
  rho=rad.mayor;
  
  mu.q= cos(theta1)*rho;
  mu.p= sin(theta1)*cos(theta2)*rho;
  xhi.q=sin(theta1)*sin(theta2)*cos(theta3)*rho;
  xhi.p=sin(theta1)*sin(theta2)*sin(theta3)*rho;
  
  for(int k=0; k<muestreo; k++){
    //Aqui hacemos la transformada de Fourier
	    // Que nos da la Funcion de Weyl
    
    weylreal+=cos((x[k].simplecticproduct(mu)+
		   y[k].simplecticproduct(xhi))/hbarra);
    
    weylimag+=sin((x[k].simplecticproduct(mu)+
		   y[k].simplecticproduct(xhi))/hbarra);
    
  };
  
  
  
  weylreal=weylreal/(hbarra*pi*muestreo);
  weylimag=weylimag/(hbarra*pi*muestreo);
	  

   
  bigrhoWeylSign=signum(weylreal);
  
  //**************************************************	  
  
  result=(bigrhoWeylSign!=smallrhoWeylSign);
  return result;
	  

};
