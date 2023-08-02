#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include "mur.hpp"

using namespace std;


complex <double> t_m(int i,complex <double> gamma_p, complex <double> gamma, double beta, double epsr, vector <double> d, double l, vector <murs> array) {


vector <int> normal = array[i].get_vecnormal(); 
double res = abs((d[0]*normal[0] + d[1]*normal[1])/(sqrt(d[0]*d[0]+d[1]*d[1]))) ;
double s = l/sqrt(1-((1/epsr)*(1-res*res))); //a changer ?
//cout << "resss est " << res << endl ; 

complex <double> exposant(0, 2*beta*s*sqrt(1.0/epsr)*(1.0-(res*res)));
complex <double> tm1 = (1.0-(gamma_p*gamma_p))*exp(-gamma*s);
complex <double> tm2 = 1.0 - (gamma_p*gamma_p)*exp(-2.0*gamma*s)*exp(exposant);
complex <double> tm = tm1/tm2;
//cout << "TM est" << tm << endl ; 
return tm ;
}
