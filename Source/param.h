/* param.h 030390 */

#include "constant.h"

#ifdef  ANSIPROTO
void SetArguments(int argc, char **argv);
BOOLEAN ExistParameter(char ch);
char    *TrueParamStr(int Item);
BOOLEAN HelpRequested(void);
int     DetermineOutputPrecision(void);
char    *NameParameter(char ch);
int     IntegerParameter(char ch, int Default, int MinVal, int MaxVal);
int     IntegerPosParam(int Item, int Default);
void	RealParameter(char ch, double *Default, double MinVal, double MaxVal);
void	RealPosParam(int Item, double *Default);
char *DisplayReal(double value, int precision,char *buffer);
BOOLEAN NullName(char *Name);
FILE *OpenInput(char *Name);
FILE *OpenOutput(char *Name);
int SizeParameter(char ch,int i,int imax);
#else
void SetArguments();
BOOLEAN ExistParameter();
char    *TrueParamStr();
BOOLEAN HelpRequested();
int     DetermineOutputPrecision();
char    *NameParameter();
int     IntegerParameter();
int     IntegerPosParam();
void    RealParameter();
void    RealPosParam();
char  *DisplayReal();
BOOLEAN NullName();
FILE *OpenInput();
FILE *OpenOutput();
int SizeParameter();
#endif

