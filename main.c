#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include "norme.hpp"
#include "mur.hpp"
#include "vecteur.hpp"
#include "coefficient_reflexion.hpp"
#include "gamma_m.hpp"
#include "t_m.hpp"
#include "comp_directe.hpp"


using namespace std;


int main()
{
    /*definition des constantes*/
	double epsr = 4.8 ;
    double pi = 3.14159265358979323846264338327950288419716939937510582;
    double f = 868.3e6 ;  // à changer en fonction du problème 
    int c = 299792458;
    double lamda = c/f;
    static double omega = 2*pi*f;
	double mu_0 = 4.0e-7*pi ; 
	double l = 0.15 ;
	
	
	vector <double> mur = {0.0,0.0}, origine = {32.0,10.0} ; 
	vector <int> directeur = {0,1} ,vec_normal = {1,0} ; 
	double beta = (2*pi)/lamda; 
	
	/*il y a trois type de murs*/

    vector <murs> array;

    for (int i = 0; i < 3; ++i) {
        array.push_back(murs(i,0.018,4.8));
    }
	
	
	complex <double> z_m = array[1].calcul_Zm(mu_0,omega) ;
    complex <double> gamma = array[1].calcul_gamma(omega, mu_0);
    
      /*PAS DE REFLEXION*/
    /*calcul de la composante directe */
    vector <double> antenne = {32.0,10.0} ;
    vector <double> point_cible = {47.0, 65.0} ; //ITERATION 
    
   
    //calcul_comp_directe(array, antenne, point_cible, z_m, gamma, epsr, beta) ;
    
    ofstream mFlux("vecteurs.dat") ;
   
     
    /*UNE REFLEXION*/
     /*calcul de gamma_perp pour la transmission et reflexion*/
	    for (int i = 0 ; i<3 ; i++) {
		//cout << "-----------------------CALCUL NEW LIGNE----------------------------------"<< endl <<endl;
	    vector <double> ri = array[i].caclul_ri(i,antenne, point_cible,array) ; 
	    vector <double> d = {point_cible[0] - ri[0] ,point_cible[1] - ri[1]};
        //cout << " RI est (" << ri[0] <<","<<ri[1] << ")" << endl ;
	    vector <double> pr = intersection(i, d, ri, array) ; // MUR DE LA REFLEXION
	    if (pr[0] != 0 || pr[1] != 0) {	
	    mFlux << " " << origine[0] << " " << origine[1] << " " << pr[0]-origine[0] << " " << pr[1]-origine[1] << endl ;	
        int numero_mur_ref = i ; 
        vector <double> vec_pr = {pr[0] - origine[0] , pr[1] - origine[1]} ; //VEC A DESSINER
        complex <double> gamma_perp5 = reflexion(i, vec_pr, z_m, epsr, array);
        complex <double> gamma_M1 = gamma_m(numero_mur_ref, gamma_perp5, gamma, d, l, epsr, beta, array) ; 	 
	       for (int j = 0 ; j<3 ; j++){
	         if (j != numero_mur_ref){
	         vector <double> vec_i = intersection(j, vec_pr, origine, array);
	             if (vec_i[0] != 0 || vec_i[1] != 0){
	             mFlux << " " << pr[0] << " " << pr[1] << " " << point_cible[0]-pr[0] << " " << point_cible[1]-pr[1]<< endl ;
	             int numero_murt = j ;
                 complex <double> gamma_perp6 = reflexion(j, vec_pr, z_m, epsr, array);
                 t_m(numero_murt, gamma_perp6, gamma, beta, epsr,vec_pr,l, array); 
                 //cout << "RESUME : on a une transmission avec le mur numéro " << numero_murt << " et une reflexion qur le mur " << numero_mur_ref << endl
                  ;} 
			 }
		 }
	 }
 }

   
    
    
    
    /*ESSAI CLASSE*/
    ofstream monFlux("murs.dat");
    for (int i = 0 ; i<3 ; i++) {
		int mat[3][2] = {{0,0}, {0,20}, {0,80}} ; 
		int mat2[3][2] = {{0,80}, {60,0}, {60,0}} ;
		array[i].sedessiner({mat2[i][0],mat2[i][1]},{mat[i][0],mat[i][1]}) ;
		}
   
    //system("gnuplot") ;
    return 0;
}
