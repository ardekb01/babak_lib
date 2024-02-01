///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
// Copyright (C) 2024 Babak A. Ardekani, PhD - All Rights Reserved.  //
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

#include <string.h>
#include <ctype.h>

int isOrientationCodeValid(const char *orientCode)
{
   char code[4];

   if( strlen(orientCode) != 3)
   {
      return(0);
   }

   strcpy(code, orientCode);

   for(int i=0; i<3; i++)
      code[i] = toupper(code[i]);

   if( 
   strcmp(code, "PIL") == 0 ||
   strcmp(code, "PIR") == 0 ||
   strcmp(code, "PSL") == 0 ||
   strcmp(code, "PSR") == 0 ||
   strcmp(code, "PLI") == 0 ||
   strcmp(code, "PLS") == 0 ||
   strcmp(code, "PRI") == 0 ||
   strcmp(code, "PRS") == 0 ||
   strcmp(code, "AIL") == 0 ||
   strcmp(code, "AIR") == 0 ||
   strcmp(code, "ASL") == 0 ||
   strcmp(code, "ASR") == 0 ||
   strcmp(code, "ALI") == 0 ||
   strcmp(code, "ALS") == 0 ||
   strcmp(code, "ARI") == 0 ||
   strcmp(code, "ARS") == 0 ||
   strcmp(code, "IPL") == 0 ||
   strcmp(code, "IPR") == 0 ||
   strcmp(code, "IAL") == 0 ||
   strcmp(code, "IAR") == 0 ||
   strcmp(code, "ILP") == 0 ||
   strcmp(code, "ILA") == 0 ||
   strcmp(code, "IRP") == 0 ||
   strcmp(code, "IRA") == 0 ||
   strcmp(code, "SPL") == 0 ||
   strcmp(code, "SPR") == 0 ||
   strcmp(code, "SAL") == 0 ||
   strcmp(code, "SAR") == 0 ||
   strcmp(code, "SLP") == 0 ||
   strcmp(code, "SLA") == 0 ||
   strcmp(code, "SRP") == 0 ||
   strcmp(code, "SRA") == 0 ||
   strcmp(code, "LPI") == 0 ||
   strcmp(code, "LPS") == 0 ||
   strcmp(code, "LAI") == 0 ||
   strcmp(code, "LAS") == 0 ||
   strcmp(code, "LIP") == 0 ||
   strcmp(code, "LIA") == 0 ||
   strcmp(code, "LSP") == 0 ||
   strcmp(code, "LSA") == 0 ||
   strcmp(code, "RPI") == 0 ||
   strcmp(code, "RPS") == 0 ||
   strcmp(code, "RAI") == 0 ||
   strcmp(code, "RAS") == 0 ||
   strcmp(code, "RIP") == 0 ||
   strcmp(code, "RIA") == 0 ||
   strcmp(code, "RSP") == 0 ||
   strcmp(code, "RSA") == 0) 
      return(1);

   return(0);
}
