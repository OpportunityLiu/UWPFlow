#pragma once

/* param.h 030390 */

#include "constant.h"

#ifdef  ANSIPROTO
void SetArguments(int argc, const char **argv);
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
BOOLEAN NullName(const char *Name);
FILE *OpenInput(const char *Name);
FILE *OpenOutput(const char *Name);
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

