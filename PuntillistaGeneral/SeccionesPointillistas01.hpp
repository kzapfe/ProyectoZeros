/*Las secciones en centros y cuerdas */

void WeylSection(vector<simplectic>& Xint, vector<simplectic>& Yint, 
		 int  maximumpoints, 
		 bool selectx, bool selecty, string nombrearchivo){
  

  
  ofstream Archivo;
  //Funcion chingona que convierte un string a c_str que funciona para open.
  Archivo.open(nombrearchivo.c_str()); 
  Archivo<<"#Sphere" <<endl;
  
  //mu es la SFT de x, xhi de y
  /*Tabla Booleana */
  //No te interesan dos cortes: los cruzados
  /* 0 0 -> muq xhiq 
     0 1 -> muq mup
     1 0 -> xhiq xhip
     1 1 -> mup xhip */
    
   double auxmuq=0.00, auxmup=0.00, auxxhiq=0.00, auxxhip=0.00;
   double coordenadauno, coordenadados;
   simplectic mu(auxmuq, auxmup), xhi(auxxhiq,auxxhip);
  
   clock_t empieza,final;
   double tiemponecesario=0.0;
   
    if((!selectx)&&(!selecty))
	    { //Corte en Q	 
	      Archivo<<"# Corte en muq, xhiq"<<endl;
	    }
	  else if((selectx)&&(!selecty))
	    { //Corte en X	      
	      Archivo<<"# Corte en muq, mup"<<endl;
	    }
	  else if ((!selectx)&&(selecty))
	    { // Corte en Y
	      Archivo<<"# Corte en xhiq,xhip"<<endl;
	    }
	  else if((selectx)&&(selecty))
	    { //Corte en P     
	      Archivo<<"# Corte en mup, xhip"<<endl;
	    }
	 


   empieza=clock();
    
      for(int n=-resol; n<resol; n++){
	for(int m=-resol; m<resol; m++){
	  //partes chatas
	  
	  double WeylReal;	
	  double WeylImag;

	  //El corte de Weyl tiene la escala AL REVEZ , aprox
	  coordenadauno=(1./magnituddeinteres)*(double)n/(double)resol;     
	  coordenadados=(1./magnituddeinteres)*(double)m/(double)resol;      
           
	  if((!selectx)&&(!selecty))
	    { //Corte en Q
	      mu.q=coordenadauno;
	      xhi.q=coordenadados;
	      
	    }
	  else if((selectx)&&(!selecty))
	    { //Corte en Mu
	      mu.q=coordenadauno;
	      mu.p=coordenadados;
	      
	    }
	  else if ((!selectx)&&(selecty))
	    { // Corte en Xhi
	      xhi.q=coordenadauno;
	      xhi.p=coordenadados;
	      
	    }
	  else if((selectx)&&(selecty))
	    { //Corte en P
	      mu.p=coordenadauno;
	      xhi.p=coordenadados;
	      
	    }
	  
	  //lets get to business
      
	  //primero las puntas diraquianas 
	  WeylReal=0.000;
	  WeylImag=0.000;
	 
	  for(int i=0; i<maximumpoints; i++){	
	    WeylReal+=cos(-(Xint.at(i).simplecticproduct(mu)+
			    Yint.at(i).simplecticproduct(xhi))/hbar);
	  
	    WeylImag+=sin(-(Xint.at(i).simplecticproduct(mu)+
			    Yint.at(i).simplecticproduct(xhi))/hbar);
	    
	  }
	  
	  WeylReal=WeylReal*2.0*pi/(maximumpoints*sqrt(hbar));
	  WeylImag=WeylImag*2.0*pi/(maximumpoints*sqrt(hbar));
	  
	  Archivo<<coordenadauno<<"\t"<<coordenadados<<"\t"
		 <<WeylReal<<"\t"<<WeylImag<<endl;
      
	}
	
    
	Archivo<<endl;
	
      } //Termina el dibujo (sobre el segundo ciclo de resoluciones)
  
      Archivo.close();
      final=clock();
      
      tiemponecesario=double(final-empieza)/CLOCKS_PER_SEC;
      cout<<"Me tarde "<<tiemponecesario
	  << "s en calcular el corte "<< 
	selectx << selecty<<
	" en Wigner. "<<endl;
      //Aqui termina su trabajo el hilo 0
      
      
      
      return;

};





