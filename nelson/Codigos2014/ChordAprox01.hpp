//Partes de la funcion de cuerdas hasta la 3 potencia


//Para mi convencion, esta mierda es J


  

double XiUno(vec centermedian, vec chords){
  //Acuerdate que las formulas las tienes en la
  //Convencion de Alfredo
  //Tu matrix simplectica aqui es al revez
  //Esta rutina solo funciona para dos g.d.l.
  mat simplecticJ=zeros(4,4);
  simplecticJ(0,1)=-1.00;
  simplecticJ(1,0)=1.00;
  simplecticJ(2,3)=-1.00;
  simplecticJ(3,2)=1.00;
  
  double result=0.000;
  
  result=as_scalar(trans(chords)*simplecticJ*centermedian);

  return result;

};
  

double XiDos(mat secundcumulant, vec chords){
  //Estaria bueno que descubrieras como hacer esto bien
  //This is NOT J, is other tingh
  mat simplecticW=eye(4,4);
  simplecticW(1,1)=-1.00;
  simplecticW(3,3)=-1.00;
  mat transpasar=zeros(4,4);
  transpasar(0,1)=1.000;
  transpasar(1,0)=1.000;
  transpasar(2,3)=1.000;
  transpasar(3,2)=1.000;


  double result=0.00;
  result=as_scalar(trans(chords)*transpasar*
		   simplecticW*secundcumulant*simplecticW*
		   transpasar*
		   chords);
  
  return result;
};


double XiTres(cube thirdcum, vec chords){
  double result=0;
  
  result=thirdcum(0,0,0)*chords(0)*chords(0)*chords(0)+
    thirdcum(1,1,1)*chords(1)*chords(1)*chords(1)+
    thirdcum(2,2,2)*chords(2)*chords(2)*chords(2)+
    thirdcum(3,3,3)*chords(3)*chords(3)*chords(3)+
    thirdcum(0,0,1)*chords(0)*chords(0)*chords(1)+
    thirdcum(0,1,1)*chords(0)*chords(1)*chords(1)+
    thirdcum(2,2,3)*chords(2)*chords(2)*chords(3)+
    thirdcum(3,3,2)*chords(3)*chords(3)*chords(2)+
    thirdcum(0,0,2)*chords(0)*chords(0)*chords(2)+
    thirdcum(2,2,0)*chords(2)*chords(2)*chords(0)+
    thirdcum(0,0,3)*chords(0)*chords(0)*chords(3)+
    thirdcum(2,2,1)*chords(2)*chords(2)*chords(1)+
    thirdcum(1,1,2)*chords(1)*chords(1)*chords(2)+
    thirdcum(3,3,0)*chords(3)*chords(3)*chords(0)+
    thirdcum(1,1,3)*chords(1)*chords(1)*chords(3)+
    thirdcum(3,3,1)*chords(3)*chords(3)*chords(1)+
    thirdcum(0,1,2)*chords(0)*chords(1)*chords(2)+
    thirdcum(2,3,0)*chords(2)*chords(3)*chords(0);

  return result;
};
