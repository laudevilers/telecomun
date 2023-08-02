#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include <complex>
#include "mur.hpp"
#include "vecteur.hpp"

using namespace std;


murs::murs(int numero_mur, double sigma, double epsr) : i(numero_mur), m_sigma(sigma), m_epsr(epsr), m_l(10)
{
}

void murs::sedessiner(vector <int> dxdy, vector <int> origine) const {
	
	ofstream monFlux("murs.dat",std::ios::app) ;
	
	if(monFlux)  //On teste si tout est OK
{
	
     monFlux << " " << origine[0] << " " << origine[1] << " " << dxdy[0] << " " << dxdy[1] << endl ;
     //system ("gnuplot") ;
}
else
{
    cout << "ERREUR: Impossible d'ouvrir le fichier." << endl ;
}
//cout << "YES" << endl;
} 


complex <double> murs::calcul_Zm( double mu_0, double omega ) const {
 
 double eps_0 = 8.85418782e-12 ; 
 double eps = m_epsr * eps_0;
 complex <double> den(eps,-m_sigma/omega);
 complex <double> zm = sqrt(mu_0/den);
 
 //cout << "z2 est" << zm << endl ;

 return zm ;
 
 }
 
 
 complex <double> murs::calcul_gamma(double omega ,double mu_0) const {
 
 double eps_0 = 8.85418782e-12 ; 
 double eps = m_epsr * eps_0;
 double beta = omega * sqrt((mu_0*eps)/2) * sqrt(sqrt(1+pow(m_sigma/(omega*eps),2))+1);
 double alpha = omega * sqrt((mu_0*eps)/2) * sqrt(sqrt(1+pow(m_sigma/(omega*eps),2))-1);
 complex <double> gamma(alpha,beta);
 //cout << "gamma est " << gamma << endl ; 
 return gamma ; 
 }
 
 
 
 int murs::get_longueur() {
	  
	int mat[] = {80, 60, 65};
    //cout << "la longueur du mur est " << mat[i] << endl ;
    return mat[i] ;
     }
    
vector <int> murs::get_origine() {
	  
	int mat2[3][2] = {{0,0}, {0,20}, {0,80}};
	vector <int> origine = {mat2[i][0],mat2[i][1]} ;
	//cout << "origine est (" << origine[0] << "," << origine[1] << ") " << endl;
    return origine ;
    }
    
    
vector <int> murs::get_vecdirecteur() {
	  
	int mat3[3][2] = {{0,1}, {1,0}, {1,0}};
	vector <int> direc = {mat3[i][0],mat3[i][1]} ;
	//cout << "direc est  (" << direc[0] << "," << direc[1] << ") " << endl;
    return direc ;
    }
    
    
vector <int> murs::get_vecnormal() {
	  
	int mat[3][2] = {{1,0}, {0,1}, {0,1}};
	vector <int> normal = {mat[i][0],mat[i][1]} ;
	//cout << "direc est  (" << direc[0] << "," << direc[1] << ") " << endl;
    return normal ;
    }
    
    
vector <double> murs::caclul_ri(int i, vector <double> antenne, vector <double> point_cible, vector <murs> array) const {
	 vector <int> vec_n = array[i].get_vecnormal() ; 
     vector <int> origine_mur = array[i].get_origine() ;
     double distance_mur = (antenne[0]-origine_mur[0])*vec_n[0] + (antenne[1]-origine_mur[1])*vec_n[1];
	 vector <double> ri = {(antenne[0]-(2*distance_mur)*vec_n[0]),(antenne[1]-(2*distance_mur)*vec_n[1])} ;
	 //vector <double> vec_ri_point_cible = {point_cible[0]-ri[0],point_cible[1]-ri[1]} ;
	 //vector <double> pr = intersection(i, vec_ri_point_cible, ri, array) ; 
	 ////cout << "VERIF" << ri[0] << " " << ri[1] << endl;
	 //if (pr[0] == 0 && pr[1] == 0) { ri = {0.0,0.0}
		 //;
		 //cout << "pas de ri"<< endl << endl ; }
	 return ri ;
}

 
 murs::~murs()
 {
	 }
