//Aqui se ponen las cuatro derivadas del potencial de Nelson
//H=(x.p^2+y.p^2)/2+x.q^2/20+(y.q-x.q^2/2)^2


double dH_dxq(simplectic &x, simplectic &y){
  double result;

  result=x.q/20.0-2.0*(y.q-x.q*x.q/2.0)*x.q;
  return result;

};


double dH_dyq(simplectic &x, simplectic &y){
  double result;

  result=2.0*(y.q-x.q*x.q/2.0);
  return result;

};


//Actually these are cumbersome, try not to use them.
//just for generality.
double dH_dxp(simplectic &x, simplectic &y){
  double result;
  result=x.p;
  return result;

};


double dH_dyp(simplectic &x, simplectic &y){
  double result;
  result=y.p;
  return result;

};
