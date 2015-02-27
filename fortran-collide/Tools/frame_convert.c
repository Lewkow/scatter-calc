#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define PI 3.14159
#define NUC 1.67 
#define ELEC 9.109e-4

int main(int argc, char *argv[]) {
	float xi, xf, m1, m2, E, ti, tf;
	float n1, n2, e1, e2;  
	printf("\n\tThis is a frame conversion program for going from Lab to CM values for collisions\n"); 
	printf("\n\tEnter projectile mass numbers [nucleons|electrons]: "); 
	scanf("%f%f", &n1, &e1); 
	printf("\tEnter target mass [nucleons|electrons]: "); 
	scanf("%f%f", &n2, &e2); 
	printf("\tEnter Lab collision energy: "); 
	scanf("%f", &E); 
	printf("\tEnter initial Lab angle (deg): "); 
	scanf("%f", &ti); 
	printf("\tEnter final Lab angle (deg): "); 
	scanf("%f", &tf); 
	ti = ti*(PI/180.0); 
	tf = tf*(PI/180.0); 	
	xi = tan(ti); 	
	xf = tan(tf); 
	float cmE, cmti, cmtf, lam;
	m1 = NUC*n1+ELEC*e1; 
	m2 = NUC*n2+ELEC*e2;  
	lam = m1/m2; 
	cmti = acos((-xi*xi*lam+sqrt(1.0+xi*xi-xi*xi*lam*lam))/(1.0+xi*xi)); 
  cmtf = acos((-xf*xf*lam+sqrt(1.0+xf*xf-xf*xf*lam*lam))/(1.0+xf*xf)); 
//	cmti = acos((-xi*xi*lam-sqrt(1.0+xi*xi-xi*xi*lam*lam))/(1.0+xi*xi)); 
//  cmtf = acos((-xf*xf*lam-sqrt(1.0+xf*xf-xf*xf*lam*lam))/(1.0+xf*xf)); 
	cmti = cmti*(180.0/PI); 
	cmtf = cmtf*(180.0/PI);
	ti = ti*(180.0/PI); 
	tf = tf*(180.0/PI);  
	cmE = m2/(m1+m2)*E; 
	printf("\n\tElab = %f -> Ecm = %f\n", E, cmE); 
	printf("\tThetaLab = %f -> ThetaCm = %f\n", ti, cmti); 
	printf("\tThetaLab = %f -> ThetaCm = %f\n\n", tf, cmtf); 

	return(0); 
}


