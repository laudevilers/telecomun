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
	//double epsr = 4.8 ;
    double pi = 3.14159265358979323846264338327950288419716939937510582;
    //double f = 868.3e6 ;  // à changer en fonction du problème 
    double f = 26e9 ;
    int c = 299792458;
    double lamda = c/f;
	//double l = 0.15 ; effacer
	//int N_total = 3;
	
	//Projet ray tracing 
	const vector <double> antenne = {-10.0,0.5};
	int N_brique = 8;
	int N_beton = 3;
	int N_cloison = 14; 
	const int N_total = N_brique+N_beton+N_cloison ;
	double max_x = 5, max_y = 5;
	double pas = 0.5, debut = 1 ;
	double gain_max = 21.5836;
	double delta = 0.1 ;
	int P_tx1 = 35;
	int ra = 73 ;
	
	
	/*vector <double> mur = {0.0,0.0}, origine = {32.0,10.0} ; 
	vector <int> directeur = {0,1} ,vec_normal = {1,0} ;*/ 
	double beta = (2*pi)/lamda; 
	complex <double> beta2(0.0,beta);
	
	/*création de la liste contenant tous les murs*/

    vector <murs> array;

    for (int i = 0; i < N_beton; ++i) {
        array.push_back(murs(i,0.014,5,0.5)); //murs(int numero_mur, double sigma, double epsr, double epaisseur)
    }
    for (int i = N_beton; i < (N_brique+N_beton) ; ++i) {
        array.push_back(murs(i,0.02,4.6,0.3)); //murs(int numero_mur, double sigma, double epsr, double epaisseur)
    }
	 for (int i =(N_brique+N_beton) ; i < (N_total) ; ++i) {
        array.push_back(murs(i,0.04,2.25,0.1)); //murs(int numero_mur, double sigma, double epsr, double epaisseur)
    }
    
    /*verification*/
    /*cout << "-----------VERIFICATION--------------"<< endl ;
    for (int i = 0 ; i < N_total ; i ++ ) {
		cout << "" << array[i].get_epsr() << endl ; }*/
    


    
      /*PAS DE REFLEXION*/
    /*calcul de la composante directe */
 
    //remove( "vecteurs.dat" );
    
    /*ouvrir un fichier pour écrire la matrice dedans et pouvoir la plot */
    fopen("matrice.dat", "w") ;
    FILE *fichier = fopen("matrice.dat", "w") ; 
 
    //EXERCICE
    /*vector <double> antenne = {32.0,10.0} ;
    vector <double> point_cible = {47.0, 65.0} ; //ITERATION */
    
    /* ===================================================
     *                    déroulement 
     * ===================================================
     * à chaque itération, on crée le vecteur direct appelé vec_direct reliant l'antenne et le point cible.
     * Ensuite, on regarde s'il y a une intersection entre ce vecteur et un des murs. Si c'est le cas un 
     * coefficient de transmission est calculé et une puissance en dBm également.  */
     
    //for (double i = debut ; i < max_x ; i ++) {
		//for (double j = debut ; j < max_y ; j ++ ) {
        //vector <double> point_cible = {i*pas, -j*pas} ; 
        ////cout << "--------------------POINT CONCERNE NUMERO--------------------------- (" << i << "," << -j << ")" << endl ;
        //vector <double> vec_direct = {point_cible[0]-antenne[0],point_cible[1]-antenne[1]};
        ////cout << "===============================================VECTEUR DIRECT :" << " vec direct est ( " << vec_direct[0] << "," << vec_direct[1] << ")" << endl ;
        ////mFlux << " " << antenne[0] << " " << antenne[1] << " " << i-antenne[0] << " " << -j-antenne[1] << endl ;
        //double norm_vec = norm(vec_direct);
        //double phi = acos(((j*pas)+antenne[1])/norm_vec);
        //double G_tx1 = gain_max - 12*((phi + (180*delta))/(30*pi)); /*en dbmt*/
        //double gain1 = (G_tx1+P_tx1)/10;
        //double gain = (pow(10,gain1))*10e-3; /*bonnes unités*/
        ////cout << "VERIFIFICATION phi  ======================================================================== " << phi  << endl;
        //complex <double> TM = calcul_comp_directe(array, vec_direct,antenne, N_total, beta, 10000) ;
        //complex <double> champs_direct = TM*sqrt(60*gain)*exp(-beta2*norm_vec)/norm_vec ;
        
        
        ////cout << "champs direct " << champs_direct << endl ;  
        //double norme_direct = norm(champs_direct) ; 
        //double prx = norme_direct*norme_direct*(lamda*lamda)/(8*(pi*pi)*ra);
        ////cout << "norme est " << norme_direct << endl ;
        //double prx_dbm = 30 + 10*log10(prx); 
        ////cout << "T_M est " << TM << endl ; 
        ////cout << "======= la puissance du point (" << i*pas << "," <<  j*pas << ")" << "est " << prx_dbm << endl ;
        //fprintf(fichier,"%lf ",prx_dbm); 
        ////complex <double> champs_direct = TM*sqrt(60*1.64e-3)*exp(-beta2*norm_vec)/norm_vec ;//EXERCICE */
        
        //}
  //fprintf(fichier,"\n") ;     
        
//}

    ///*fermeture du fichier*/
    //fclose(fichier) ; 

      
    /* ===================================================
     *                    déroulement 
     * ===================================================
     * à chaque itération, on calcule l'image de l'antenne avec le mur i. On crée le vecteur d reliant 
     * le point_cible (iteration sur les points) et ri. On regarde ensuite s'il y a une intersection de ce vecteur avec l'un d'un mur
     * (iteration sur tous les murs) .Si c'est le cas, il y aura une reflexion possible sur ce mur. On calcule alors l'intersection de 
     * d avec le mur, c'est notre point pr. Il faut maintenant trouver le mur avec lequel il y a transmission. on cree le vecteur vec_pr 
     * et on cherche avec quel mur il y a une transmission. pour avoir le coeffiient total il faut multiplier le coefficient de transmissin
     * et de reflexion  */   
     
     
     
    /*UNE REFLEXION*/
     /*calcul de gamma_perp pour la transmission et reflexion*/
	    for (int i = debut ; i < max_x ; i++) {
			for (int j = debut ; j< max_y ; j++) {		  
		cout << "-----------------------CALCUL NEW LIGNE----------------------------------"<< endl <<endl;
		vector <double> point_cible = {i*pas, -j*pas} ; 
		vector <double> vec_direct = {point_cible[0]-antenne[0],point_cible[1]-antenne[1]};
		
		for (int mur = 0 ; mur < N_total ; mur ++){
		/*CALCUL RI*/	
	    vector <double> ri = array[mur].caclul_ri(mur,antenne, point_cible,array) ; 
	    vector <double> d = {point_cible[0] - ri[0] ,point_cible[1] - ri[1]};
	    vector <double> pr = intersection(mur, d, ri, array) ; // MUR DE LA REFLEXION //intersection de d avec le mur i ?
	    /*SI ON A UNE INTERSECTION ENTRE RI ET POINT CIBLE = REFLEXION */
	    if (pr[0] != 0 || pr[1] != 0) {	
		cout << " RI est (" << ri[0] <<","<<ri[1] << ")" << endl ;	
	    //mFlux << " " << origine[0] << " " << origine[1] << " " << pr[0]-origine[0] << " " << pr[1]-origine[1] << endl ;//EXERCICE
        int numero_mur_ref = mur ; 
        vector <double> vec_pr = {pr[0] - antenne[0] , pr[1] - antenne[1]} ; //VEC A DESSINER
        complex <double> gamma_perp5 = reflexion(mur, vec_pr, array); 
        complex <double> gamma_M1 = gamma_m(numero_mur_ref, gamma_perp5, vec_pr, beta, array) ; 	
        complex <double> TM = calcul_comp_directe(array, vec_pr ,antenne, N_total, beta, numero_mur_ref) ; 
        //cout << " TM est " << TM << " et GAMMA M est " << gamma_M1 << endl ;
        double norm_vec = norm(d); 
	    double norme2 = norm(vec_direct);
        double phi = acos(((j*pas)+antenne[1])/norme2); //pas norme vec ici c norme vec direct
        double G_tx1 = gain_max - 12*((phi + (180*delta))/(30*pi)); /*en dbmt*/
        double gain1 = (G_tx1+P_tx1)/10;
        double gain = (pow(10,gain1))*10e-3; /*bonnes unités*/
        // A T ON LE M1 et TM final?
	    complex <double> champs = gamma_M1*TM*sqrt(60*gain)*exp(-beta2*norm_vec)/norm_vec ;
	    complex <double> champs_total = champs_total +  champs; //Somme des champs  
	     cout <<"champs est " << champs << " et champs total est " << champs_total << endl ;     

        //fprintf(fichier,"%lf ",prx2_dbm); 
	    
	       /*for (int j = 0 ; j<3 ; j++){
	         if (j != numero_mur_ref){
	         vector <double> vec_i = intersection(j, vec_pr, origine, array);
	             if (vec_i[0] != 0 || vec_i[1] != 0){
	             mFlux << " " << pr[0] << " " << pr[1] << " " << point_cible[0]-pr[0] << " " << point_cible[1]-pr[1]<< endl ; //EXERCICE
	             int numero_murt = j ;
                 complex <double> gamma_perp6 = reflexion(j, vec_pr, z_m, epsr, array);
                 t_m(numero_murt, gamma_perp6, gamma, beta, epsr,vec_pr,l, array); 
                 cout << "RESUME : on a une transmission avec le mur numéro " << numero_murt << " et une reflexion qur le mur " << numero_mur_ref << endl
                  ;} 
			 }
		 }*/
	 }
	         
   }
        cout <<" !!!!! champs est " << champs_total << " et champs total est " << champs_total << endl ; 
        double norme_champs = norm(champs_total) ; 
        double prx2 = norme_champs*norme_champs*(lamda*lamda)/(8*(pi*pi)*ra);
        cout << "norme est " << norme_champs << endl ;
        double prx2_dbm = 30 + 10*log10(prx2); 
        cout << "======= la puissance du point (" << i*pas << "," <<  j*pas << ")" << "est " << prx2_dbm << endl ;
 }
 //fprintf(fichier,"\n") ; 
}

    /*fermeture du fichier*/
    fclose(fichier) ;

   
    
    /*ESSAI CLASSE*/
    ofstream monFlux("murs.dat");
    for (int i = 0 ; i < N_total ; i++) {
		//int mat[N_total][2] = {{0,0}, {0,20}, {0,80}} ; 
		//int mat2[N_total][2] = {{0,80}, {60,0}, {60,0}} ;
		
		double mat[N_total][2] = {{0,-45}, {75,-45}, {75,-70},{0,0}, {0,0}, {100,0},{85,0}, {85,-14.8398}, {85,-27.6795},{93.5,-27.6795}, {85,-27.6795}, {15,0},
		{0,-9}, {15,-5}, {0,-18},{15,-14}, {0,-27}, {15,-23},{0,-36}, {15,-32}, {15,-41},{40,-15}, {50,-15}, {35,-30}, {55,-30}};
		
		vector <double> dxdy = array[i].get_dxdy(i) ;
			
		array[i].sedessiner(dxdy,{mat[i][0],mat[i][1]}) ; //mettre un get_origine?
		}
   
    system("gnuplot -persistent 'f.txt'") ;
    return 0;
}
