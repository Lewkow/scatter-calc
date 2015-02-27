#include "systemheaders.h"

int main() {
  int type, in, fn;
  float val, valpow; 
  printf("\n\tUnit Conversion Program\n");
  printf("\n\tPlease select type of unit to convert:\n"); 
  printf("\t\t1 => Length\n"); 
  printf("\t\t2 => Mass\n"); 
  printf("\t\t3 => Energy\n"); 
  scanf("%d", &type); 
  if (type == 1) {
      double CMtoM, AtoM, BOHRtoM; 
      CMtoM = 1.0e-2; 
      AtoM = 1.0e-10; 
      BOHRtoM = 52.9177e-12; 
      
      printf("\tLength Conversion Selected\n");
      printf("\n\tPlease select original units:\n"); 
      printf("\t\t1 => Meter [m]\n"); 
      printf("\t\t2 => Centimeter [cm]\n"); 
      printf("\t\t3 => Angstrom [A]\n"); 
      printf("\t\t4 => Bohr Radius [a0]\n"); 
      scanf("%d", &in);
      printf("\t%d selected\n", in); 

      printf("\n\tPlease select final units:\n"); 
      printf("\t\t1 => Meter [m]\n"); 
      printf("\t\t2 => Centimeter [cm]\n"); 
      printf("\t\t3 => Angstrom [A]\n"); 
      printf("\t\t4 => Bohr Radius [a0]\n"); 
      scanf("%d", &fn);
      printf("\t%d selected\n", fn); 

      printf("\n\tPlease enter value to convert:\n"); 
      scanf("%f", &val);
      printf("\tInitial Value: %f\n", val); 

      printf("\n\tPlease enter power to convert:\n"); 
      scanf("%f", &valpow); 
      printf("\tPower to convert: %f\n", valpow); 

      if (in == 2) {
          val = val*pow(CMtoM,valpow); 
      }
      else if (in == 3) {
          val = val*pow(AtoM,valpow); 
      }
      else if (in == 4) {
          val = val*pow(BOHRtoM,valpow); 
      }
      else printf("You have input an invalid original unit\n"); 

      if (fn == 2) val = val/pow(CMtoM,valpow); 
      else if (fn == 3) val = val/pow(AtoM,valpow); 
      else if (fn == 4) val = val/pow(BOHRtoM,valpow); 

      printf("\n\tConverted Value: %e\n", val); 
  
  }

  return 0; 

}
